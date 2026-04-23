#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <assert.h>
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop
#include "yak-priv.h"
#include "ketopt.h"
#include "bseq.h"
#include <map>
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include "kseq.h" // FASTA/Q parser
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <inttypes.h>
using namespace std;

KSEQ_INIT(gzFile, gzread)
#define CHUNK_SIZE 20000000
#define min_num_shapmers 3
std::unordered_map<int, int> node_type_map1;
std::unordered_map<uint64_t, int> g_read_length_map1;
std::unordered_map<uint64_t, int> g_read_length_map_hap1;
std::unordered_map<uint64_t, int> g_read_length_map_hap2;

vector<char *> h1tg_names;
vector<char *> h2tg_names;
typedef struct
{
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
	uint64_t *r;
	bool *f;
	uint32_t *pos; // 新增：记录kmer首个碱基在序列中的位置
} ch_buf_t;

void printBits(size_t const size, void const *const ptr)
{
	unsigned char *b = (unsigned char *)ptr;
	unsigned char byte;
	int i, j;

	for (i = size - 1; i >= 0; i--)
	{
		for (j = 7; j >= 0; j--)
		{
			byte = (b[i] >> j) & 1;
			printf("%u", byte);
		}
	}
	puts("");
}

// 一个桶 b 其实是一个 动态数组，数组大小 = b->m，当前已用 = b->n。

// 桶里每条记录存了一行信息：

// a：k-mer 的哈希值（uint64_t）

// f：这个 k-mer 是 forward 还是 reverse（布尔）

// pos：k-mer 在原序列中的位置（uint32_t）

// r：这个 k-mer 属于哪条 read（seq_id）
static inline void ch_insert_buf(ch_buf_t *buf, int p, uint64_t y, uint64_t seq_id, bool f, uint32_t kmer_pos)
{
	int pre = y & ((1 << p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m)
	{
		b->m = b->m < 8 ? 8 : b->m + (b->m >> 1);
		REALLOC(b->a, b->m);
		REALLOC(b->r, b->m);
		REALLOC(b->f, b->m);
		REALLOC(b->pos, b->m); // 新增
	}
	b->a[b->n] = y;
	b->f[b->n] = f;

	b->pos[b->n] = kmer_pos; // 新增
	b->r[b->n++] = seq_id;
}

static void count_seq_buf(ch_buf_t *buf, int k, int p, int len, const char *seq, uint32_t seq_id) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	bool f;
	uint64_t x[2], mask = (1ULL << k * 2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i)
	{
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4)
		{												   // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;				   // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift; // reverse strand
			if (++l >= k)
			{ // we find a k-mer
				f = x[0] < x[1];
				uint64_t y = f ? x[0] : x[1];
				ch_insert_buf(buf, p, yak_hash64(y, mask), seq_id, f, i - k + 1); // 新增参数
			}
		}
		else
			l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

static void count_seq_buf_long(ch_buf_t *buf, int k, int p, int len, const char *seq, uint16_t seq_id) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	bool f;
	uint64_t x[4], mask = (1ULL << k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i)
	{
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4)
		{ // not an "N" base
			x[0] = (x[0] << 1 | (c & 1)) & mask;
			x[1] = (x[1] << 1 | (c >> 1)) & mask;
			x[2] = x[2] >> 1 | (uint64_t)(1 - (c & 1)) << shift;
			x[3] = x[3] >> 1 | (uint64_t)(1 - (c >> 1)) << shift;
			uint64_t y = f ? x[0] : x[1];
			if (++l >= k)
				ch_insert_buf(buf, p, yak_hash64(y, mask), seq_id, f, i - k + 1); // 新增参数
		}
		else
			l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
}

typedef struct
{ // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	int create_new;
	kseq_t *ks;
	yak_ch_t *h;
	yak_ch_t *h_pos; // 新增：专门存pos的哈希表
	uint64_t global_counter;
	int seg_n;
	std::vector<char *> *names;
	std::unordered_map<uint64_t, int> g_read_length_map;
} pldat_t;

typedef struct
{ // data structure for each step in kt_pipeline()
	pldat_t *p;
	int n, m, sum_len, nk;
	uint64_t global_bias;
	int *len;
	char **seq;
	char **qual;
	ch_buf_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
	stepdat_t *s = (stepdat_t *)data;
	yak_ch_t *h = s->p->h;
	yak_ch_t *h_pos = s->p->h_pos; // 新增
	ch_buf_t *b = &s->buf[i];
	// b->n_ins += yak_ch_insert_list_kmer_record_mapping(h, s->p->create_new, b->n, b->a, b->f, NULL, b->r, i, NULL);
	//    b->n_ins += yak_ch_insert_list_kmer_record_mapping(
	//        h, s->p->create_new, b->n, b->a, b->f, b->pos, b->r, i, NULL
	//    );
	//    插入pos表
	//    cout<<"in yak_ch_insert_list_kmer_pos"<<endl;
	//    b->n_ins += yak_ch_insert_list_kmer_pos(h, h_pos, s->p->create_new, b->n, b->a, b->r, b->pos, i);
	b->n_ins += yak_ch_insert_list_kmer_full(h, h_pos, s->p->create_new, b->n, b->a, b->r, b->pos, b->f, i);
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
	pldat_t *p = (pldat_t *)data;
	if (step == 0)
	{ // step 1: read a block of sequences
		int ret;
		stepdat_t *s;
		CALLOC(s, 1);
		s->p = p;
		s->global_bias = p->global_counter;
		while ((ret = kseq_read(p->ks)) >= 0)
		{
			p->global_counter++;
			int l = p->ks->seq.l;
			if (l < p->opt->k)
				continue;
			if (s->n == s->m)
			{
				s->m = s->m < 16 ? 16 : s->m + (s->n >> 1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);

			p->names->push_back(strdup(p->ks->name.s));

			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->opt->k + 1;
			if (s->sum_len >= p->opt->chunk_size)
				break;
		}
		if (s->sum_len == 0)
			free(s);
		else
			return s;
	}
	else if (step == 1)
	{ // step 2: extract k-mers
		stepdat_t *s = (stepdat_t *)in;
		int i, n = 1 << p->opt->pre, m;
		CALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i)
		{
			s->buf[i].m = m;
			MALLOC(s->buf[i].a, m);
			MALLOC(s->buf[i].r, m);
			MALLOC(s->buf[i].f, m);
			MALLOC(s->buf[i].pos, m); // 新增：分配pos数组
		}
		for (i = 0; i < s->n; ++i)
		{
			if (p->opt->k < 32)
				count_seq_buf(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i], s->global_bias + i);
			else
				count_seq_buf_long(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i], s->global_bias + i);
			free(s->seq[i]);
		}
		free(s->seq);
		free(s->len);
		return s;
	}
	else if (step == 2)
	{ // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t *)in;
		int i, n = 1 << p->opt->pre;
		uint64_t n_ins = 0;

		kt_for(p->opt->n_thread, worker_for, s, n);
		for (i = 0; i < n; i++)
		{
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
			free(s->buf[i].r);
			free(s->buf[i].f);
			free(s->buf[i].pos); // 新增：释放pos数组
		}
		p->h->tot += n_ins;
		free(s->buf);
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences; %ld distinct k-mers in the hash table\n", __func__,
				yak_realtime(), yak_cputime() / yak_realtime(), s->n, (long)p->h->tot);
		fprintf(stderr, "Total Processed %ld sequences\n", p->global_counter);
		free(s);
	}
	return 0;
}

pldat_t *yak_count_multi_new(const char *fn, const yak_copt_t *opt, yak_ch_t *h0)
{
	pldat_t *pl;
	CALLOC(pl, 1);
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0)
		return 0;
	pl->ks = kseq_init(fp);
	pl->opt = opt;
	if (h0)
	{
		pl->h = h0, pl->create_new = 0;
		assert(h0->k == opt->k && h0->pre == opt->pre);
	}
	else
	{
		pl->create_new = 1;
		pl->h = yak_ch_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	}
	pl->h_pos = yak_ch_init(opt->k, 30, opt->bf_n_hash, opt->bf_shift); // 新增：初始化存pos的哈希表
	pl->names = new std::vector<char *>();

	kt_pipeline(3, worker_pipeline, pl, 3);
	kseq_destroy(pl->ks);
	gzclose(fp);

	return pl;
}

// Map unitigs to hic data.

typedef struct
{
	int c[16];
	int sc[2];
	int nk;
} tb_cnt_t;

typedef struct
{
	int max;
	uint32_t *s;
} tb_buf_t;

typedef struct
{
	// char *name;
	uint16_t unit_id;
	bool forward;
	// uint64_t pos;
} mapping_res_t;

typedef struct
{
	int k, n_threads, print_diff;
	double ratio_thres;
	bseq_file_t *fp;
	yak_ch_t *ch;
	yak_ch_t *h_pos; // 新增
	tb_buf_t *buf;
	std::vector<mapping_res_t> *mappings;
	int record_num;
} tb_shared_t;

typedef struct
{
	int k, n_threads, print_diff;
	double ratio_thres;
	bseq_file_t *fp;
	yak_ch_t *ch1;
	yak_ch_t *ch2;
	yak_ch_t *h_pos1; // 新增
	yak_ch_t *h_pos2; // 新增
	tb_buf_t *buf;
	std::vector<mapping_res_t> *mappings;
	vector<char *> *contigs_names_hap1;
	vector<char *> *contigs_names_hap2;
	int record_num;
} tb_shared_hifi_t;

typedef struct
{
	int n_seq;
	tb_shared_t *aux;
	bseq1_t *seq;
	uint16_t *mappings;
	bool *mappings_forward;
	uint32_t *map_pos;
} tb_step_t;

typedef struct
{
	int n_seq;
	tb_shared_hifi_t *aux;
	bseq1_t *seq;
	uint16_t *mappings;
	bool *mappings_forward;
	uint32_t *map_pos;
} tb_step_hifi;

std::string decode_kmer(uint64_t x, int k)
{
    static const char nt4[] = {'A', 'C', 'G', 'T'};
    std::string kmer(k, 'A'); // 初始化为长度为 k 的字符串
    for (int i = 0; i < k; ++i)
    {
        kmer[k - i - 1] = nt4[x & 0x3]; // 取出最后两位并转换为字符
        x >>= 2;                        // 向右移动两位，处理下一对二进制位
    }
    return kmer;
}

struct ReadState
{
	int last_q_pos = -1; // 上一次 query kmer 的 start pos（在 query read 中）
	int last_t_pos = -1; // 上一次 target pos（在目标 read 中）
	bool has_last = false;
	bool last_forward = false;
	int consist_count = 0; // 本 read 上被判为一致（monotonic + distance ok）的次数
};

typedef struct
{
	std::string kmer;
	int qposi; // position
	int tposi;
	int forward;
	uint16_t target_id;

} Hamper_kmer;

typedef struct
{
	// std::string hapmer;
	uint64_t hapmer;
	int qposi;
	int tposi;
	int strand;

} Tread2shapmers;
static void tb_worker(void *_data, long k, int tid)
{
	tb_step_t *t = (tb_step_t *)_data;
	tb_shared_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	tb_buf_t *b = &aux->buf[tid];
	uint64_t x[4], mask;
	int i, l, shift;
	bool is_forward;
	std::map<uint16_t, int> counts;
	std::map<uint16_t, bool> forward;
	std::map<uint16_t, uint32_t> positions;
	yak_ch_t *h_pos = aux->h_pos; // 新增
	// cout << "Processing read: " << s->name << endl;
	std::map<int, std::vector<Tread2shapmers>> read2hapmers;
	if (aux->ch->k < 32)
	{
		mask = (1ULL << 2 * aux->ch->k) - 1;
		shift = 2 * (aux->ch->k - 1);
	}
	else
	{
		mask = (1ULL << aux->ch->k) - 1;
		shift = aux->ch->k - 1;
	}
	if (s->l_seq > b->max)
	{
		b->max = s->l_seq;
		kroundup32(b->max);
		b->s = (uint32_t *)realloc(b->s, b->max * sizeof(uint32_t));
	}
	memset(b->s, 0, s->l_seq * sizeof(uint32_t));
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->l_seq; ++i)
	{
		int c = seq_nt4_table[(uint8_t)s->seq[i]];
		uint32_t res;
		if (c < 4)
		{
			if (aux->ch->k < 32)
			{
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			}
			else
			{
				x[0] = (x[0] << 1 | (c & 1)) & mask;
				x[1] = (x[1] << 1 | (c >> 1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c & 1)) << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c >> 1)) << shift;
			}
			if (++l >= aux->k)
			{
				int type = 0, c1, c2;
				uint64_t y, y2;
				// cout << "kmer=" << decode_kmer(x[0], aux->ch->k) << endl;
				if (aux->ch->k < 32)
				{
					is_forward = x[0] < x[1];
					y = yak_hash64(is_forward ? x[0] : x[1], mask);
					y2 = yak_hash64(is_forward ? x[1] : x[0], mask);
				}
				else
				{
					y = yak_hash_long(x);
				}
				res = yak_ch_get_k(aux->ch, y);
				uint16_t utg = res & YAK_KEY_MASK;
				uint64_t kmer = res >> YAK_COUNTER_BITS;
				if (res != -1)
				{
					// printf("0x%016" PRIx64 "\n", y);
					//  cout << "kmer inside=" << decode_kmer(x[0], aux->ch->k) << std::endl;
					//  cout << "backward kmer inside=" << decode_kmer(x[1], aux->ch->k) << std::endl;
					// cout << "utg=" << utg << endl;
					uint32_t pos = -1;
					// cout<<y<<endl;
					uint16_t read_id = yak_ch_get_pos(aux->ch, aux->h_pos, y, &pos);
					// cout<<"is_forward: "<<is_forward<<", stored_forward_bit: "<<((res & YAK_FORWARD_MASK) != 0)<<endl;
					//  TODO:
					bool stored_forward_bit = (res & YAK_FORWARD_MASK) != 0;
					bool is_main_forward = !(stored_forward_bit ^ is_forward);
					int klen = aux->ch->k;
					int q_pos = i - klen + 1; // current k-mer start in query read

					// TODO:
					if (pos != -1) // pos == 65535 means no position info
					{
						int target_id = res & YAK_KEY_MASK; // target_id

						Tread2shapmers hap;
						hap.hapmer =x[0];					  // query kmer 的序列
						hap.qposi = l - aux->ch->k;									  // query kmer 在 read 上的位置
						hap.tposi = pos;											  // target 上的位置（来自 yak_ch_get_pos）
						hap.strand = !(((res & YAK_FORWARD_MASK) != 0) ^ is_forward); // 相同为0，不同为1
						read2hapmers[target_id].push_back(hap);
						// printf("Found hapmer: %s, QPos: %d, TPos: %d, Strand: %d, TargetID: %d\n",
						// 	   hap.hapmer.c_str(), hap.qposi, hap.tposi, hap.strand, target_id);
						// printf("utg=%u, read_id=%u, forward=%d, position=%u\n",
						// 	   utg, read_id, is_forward, pos);
					}
				}
			}
		}
		else
			l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}

	// TODO:

	for (const auto &pair : read2hapmers)
	{
		int target_id = pair.first;
		const std::vector<Tread2shapmers> &shapmers = pair.second;
		int num_shapmers = shapmers.size();
		// std::cout << "Target ID: " << target_id << "\n";
		// for (const auto &h : shapmers)
		// {
		// 	std::cout << "  Hapmer: " << h.hapmer
		// 			  << ", QPos: " << h.qposi
		// 			  << ", TPos: " << h.tposi
		// 			  << ", Strand: " << h.strand << "\n";
		// }

		if (num_shapmers < min_num_shapmers)
		{
			continue;
		}
		// 计算相同 strand 和不同 strand 的数量
		int num_same_strand = 0;
		int num_diff_strand = 0;

		// 遍历 shapmers，统计相同和不同 strand 的数量
		for (const Tread2shapmers &shapmer : shapmers)
		{
			if (shapmer.strand == 1)
			{
				num_same_strand++;
			}
			else
			{
				num_diff_strand++;
			}
		}
		// 确定是同向比对，还是反向互补的比对上
		int strand = (num_same_strand > num_diff_strand) ? 1 : 0;

		std::vector<Tread2shapmers> strand_correct_shapmers;
		// 遍历 shapmers，只挑出正确比对方向的hamper，过滤掉错误噪音
		for (const Tread2shapmers &shapmer : shapmers)
		{
			if (shapmer.strand == strand)
			{
				strand_correct_shapmers.push_back(shapmer);
			}
		}
		// 统计符合 strand 的 hapmers 数量
		int num_strand_correct_shapmers = strand_correct_shapmers.size();
		if (num_strand_correct_shapmers < min_num_shapmers)
		{
			continue;
		}
		std::vector<Tread2shapmers> posi_correct_shapmers;
		int n = 3;		 // 使用 n 个下一个 hapmers 来判断当前 hapmer
		int shift = 200; // 允许的最大位置偏移
		if (num_strand_correct_shapmers >= 2 * n)
		{
			if (strand == 1) // 相同的 strand
			{
				for (int i = 0; i <= num_strand_correct_shapmers - n - 1; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i + j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) <= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
						{
							score += -1;
						}
						else
						{
							score += 1;
						}
					}

					if (score > 0) // 认为这个 hapmer 是可信的，否则移除
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}

				// 处理最后几个 hapmers（这些 hapmers 没有 n 个后续 hapmers 来判断）
				for (int i = num_strand_correct_shapmers - n; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i - j];
						int qposi_p = hapmer_j.qposi;
						int tposi_p = hapmer_j.tposi;

						if ((qposi - qposi_p) <= 0 || (tposi - tposi_p) <= 0 || std::abs(std::abs(qposi - qposi_p) - std::abs(tposi - tposi_p)) > shift)
						{
							score += -1;
						}
						else
						{
							score += 1;
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
			//
			else // 不相同的 strand
			{
				for (int i = 0; i <= num_strand_correct_shapmers - n - 1; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i + j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) >= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
						{
							score += -1;
							// cout<<"Now score -1 : "<<score<<endl;
						}
						else
						{
							score += 1;
							// cout<<"Now score +1 : "<<score<<endl;
						}
					}

					if (score > 0) // 认为这个 hapmer 是可信的，否则移除
					{
						// cout<<"Accepted hapmer: "<<strand_correct_shapmers[i].hapmer<<", QPos: "<<strand_correct_shapmers[i].qposi<<", TPos: "<<strand_correct_shapmers[i].tposi<<", Strand: "<<strand_correct_shapmers[i].strand<<endl;
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}

				// 处理最后几个 hapmers（这些 hapmers 没有 n 个后续 hapmers 来判断）
				for (int i = num_strand_correct_shapmers - n; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i - j];
						int qposi_p = hapmer_j.qposi;
						int tposi_p = hapmer_j.tposi;

						if ((qposi - qposi_p) <= 0 || (tposi - tposi_p) >= 0 || std::abs(std::abs(qposi - qposi_p) - std::abs(tposi - tposi_p)) > shift)
						{
							score += -1;
						}
						else
						{
							score += 1;
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
		}
		else if (num_strand_correct_shapmers >= 2 && num_strand_correct_shapmers < 2 * n)
		{
			if (strand == 1) // 相同的 strand
			{
				for (int i = 0; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 0; j < num_strand_correct_shapmers; ++j)
					{
						if (i == j)
						{
							continue;
						}

						const auto &hapmer_j = strand_correct_shapmers[j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if (j > i)
						{
							if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) <= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
						else
						{
							if ((qposi_n - qposi) >= 0 || (tposi_n - tposi) >= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
			else // 不同的 strand
			{
				for (int i = 0; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 0; j < num_strand_correct_shapmers; ++j)
					{
						if (i == j)
						{
							continue;
						}

						const auto &hapmer_j = strand_correct_shapmers[j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if (j > i)
						{
							if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) >= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
						else
						{
							if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) <= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
		}
		else
		{
			continue;
		}
		if (posi_correct_shapmers.size() < min_num_shapmers)
		{
			continue;
		}
		int num_nonovlp_shapmers = 1; // 初始值为 1，考虑第一个 hapmer
		for (size_t i = 0; i < posi_correct_shapmers.size() - 1; ++i)
		{
			const Tread2shapmers &hapmer1 = posi_correct_shapmers[i];
			const Tread2shapmers &hapmer2 = posi_correct_shapmers[i + 1];

			int qposi1 = hapmer1.qposi;
			int tposi1 = hapmer1.tposi;
			int qposi2 = hapmer2.qposi;
			int tposi2 = hapmer2.tposi;
			// cout << "  Filtered Hapmer: " << hapmer1.hapmer
			// 	 << ", QPos: " << hapmer1.qposi
			// 	 << ", TPos: " << hapmer1.tposi
			// 	 << ", Strand: " << hapmer1.strand << "\n";

			// 检查相邻 hapmers 是否满足距离要求
			if (std::abs(qposi2 - qposi1) >= aux->ch->k && std::abs(tposi2 - tposi1) >= aux->ch->k)
			{
				num_nonovlp_shapmers++;
			}
		}
		if (posi_correct_shapmers.size() > 20)
		{
			int num = posi_correct_shapmers.size() / 20;
			num_nonovlp_shapmers += num;
		}
		// cout << "Num non-overlapping shapmers: " << num_nonovlp_shapmers << "\n";
		counts[target_id] = num_nonovlp_shapmers;
		forward[target_id] = (strand == 1);
		positions[target_id] = posi_correct_shapmers[0].tposi;
	}

	uint16_t max_idx = -1;
	int max_val = 0;
	for (auto idx : counts)
	{
		if (max_val < idx.second)
		{
			max_val = idx.second;
			max_idx = idx.first;
		}
	}
	t->mappings[k] = max_idx;
	t->mappings_forward[k] = forward[max_idx];
	// t->map_pos[k] = positions[max_idx];

}



static void *tb_pipeline(void *shared, int step, void *_data)
{
	tb_shared_t *aux = (tb_shared_t *)shared;
	if (step == 0)
	{
		tb_step_t *s;
		s = (tb_step_t *)calloc(1, sizeof(tb_step_t));
		s->seq = bseq_read(aux->fp, CHUNK_SIZE, 0, &s->n_seq);
		s->aux = aux;
		if (s->n_seq)
		{
			s->mappings = (uint16_t *)calloc(s->n_seq, sizeof(uint16_t));
			s->mappings_forward = (bool *)calloc(s->n_seq, sizeof(bool));
			s->map_pos = (uint32_t *)calloc(s->n_seq, sizeof(uint32_t));
			// fprintf(stderr, "[M::%s] read %d sequences\n", __func__, s->n_seq);
			return s;
		}
		else
			free(s);
	}
	else if (step == 1)
	{
		int i;
		tb_step_t *s = (tb_step_t *)_data;
		kt_for(aux->n_threads, tb_worker, s, s->n_seq);
		for (i = 0; i < s->n_seq; ++i)
		{
			mapping_res_t res;
			// res.name = strdup(s->seq[i].name);
			res.unit_id = s->mappings[i];
			res.forward = s->mappings_forward[i];
			// res.pos = s->map_pos[i];
			aux->mappings->push_back(res);
			free(s->seq[i].name);
			free(s->seq[i].seq);
			free(s->seq[i].qual);
			free(s->seq[i].comment);
		}
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences;\n", __func__,
				yak_realtime(), yak_cputime() / yak_realtime(), s->n_seq);
		free(s->seq);
		free(s->mappings);
		free(s);
	}
	return 0;
}
void do_mapping2(pldat_t *pl, bseq_file_t *hic_fn1, bseq_file_t *hic_fn2, char *out_fn, std::unordered_map<int, int> &node_type_map1)
{
	if (hic_fn1 == 0 || hic_fn2 == 0)
	{
		fprintf(stderr, "ERROR: Please give two hic files\n");
		exit(1);
	}

	ketopt_t o = KETOPT_INIT;
	int i, c, min_cnt = 2, mid_cnt = 5;
	tb_shared_t aux;

	memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = pl->opt->n_thread, aux.print_diff = 0;
	aux.ratio_thres = 0.33;
	aux.k = pl->h->k;
	aux.ch = pl->h;
	aux.h_pos = pl->h_pos; // 新增
	aux.record_num = pl->global_counter;
	cout << "Start mapping, total unitig number: " << aux.record_num << endl;
	aux.fp = hic_fn1;
	std::vector<mapping_res_t> *map1 = new std::vector<mapping_res_t>();
	aux.mappings = map1;
	aux.buf = (tb_buf_t *)calloc(aux.n_threads, sizeof(tb_buf_t));
	cout << "Start hic1 mapping " << endl;

		//***************************	 */
	long peak_bytes = yak_peakrss();
	double peak_mb = peak_bytes / 1024.0 / 1024.0;
	double peak_gb = peak_bytes / 1024.0 / 1024.0 / 1024.0;
	printf("\n[Memory Info]\n");
    printf("Peak RSS: %ld Bytes\n", peak_bytes);
    printf("Peak RSS: %.2f MB\n", peak_mb);
    printf("Peak RSS: %.2f GB\n", peak_gb);

	//********************************	 */
	kt_pipeline(2, tb_pipeline, &aux, 2);

	peak_bytes = yak_peakrss();
	peak_mb = peak_bytes / 1024.0 / 1024.0;
	peak_gb = peak_bytes / 1024.0 / 1024.0 / 1024.0;
	printf("\n[Memory Info]\n");
    printf("Peak RSS: %ld Bytes\n", peak_bytes);
    printf("Peak RSS: %.2f MB\n", peak_mb);
    printf("Peak RSS: %.2f GB\n", peak_gb);
	bseq_close(aux.fp);
	for (i = 0; i < aux.n_threads; ++i)
	{
		free(aux.buf[i].s);
	}
	free(aux.buf);

	aux.fp = hic_fn2;
	std::vector<mapping_res_t> *map2 = new std::vector<mapping_res_t>();
	aux.mappings = map2;
	aux.buf = (tb_buf_t *)calloc(aux.n_threads, sizeof(tb_buf_t));
	cout << "Start hic2 mapping " << endl;
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	for (i = 0; i < aux.n_threads; ++i)
	{
		free(aux.buf[i].s);
	}
	free(aux.buf);

	printf("Map1 size: %ld;\n", map1->size());
	printf("Map2 size: %ld;\n", map2->size());

	uint64_t success_counter_result = 0;
	uint32_t **connections_forward;
	CALLOC(connections_forward, pl->global_counter);
	for (int i = 0; i < pl->global_counter; i++)
	{
		CALLOC(connections_forward[i], pl->global_counter);
		memset(connections_forward[i], 0, sizeof(connections_forward[i]));
	}
	uint32_t **connections_backward;
	CALLOC(connections_backward, pl->global_counter);
	for (int i = 0; i < pl->global_counter; i++)
	{
		CALLOC(connections_backward[i], pl->global_counter);
		memset(connections_backward[i], 0, sizeof(connections_backward[i]));
	}
	uint32_t coverage[pl->global_counter];
	memset(coverage, 0, sizeof(uint32_t) * pl->global_counter);
	std::vector<uint64_t> failed_matches;
	const uint16_t max_count = -1;
	int map1_usefulCount = 0;
	int map2_usefulCount = 0;
	int skipped_count = 0;
	uint64_t real_useful = 0;

	//......................
	// std::ofstream outFile("output.csv");
    // if (!outFile.is_open()) {
    //     std::cerr << "Failed to open output.csv for writing.\n";
    //     return;
    // }
	// outFile << "utg_name1,pos1,utg_name2,pos2\n";
	//......................


	for (uint64_t j = 0; j < std::min(map1->size(), map2->size()); j++)
	{
		if ((*map1)[j].unit_id != 65535)
			map1_usefulCount++;
		if ((*map2)[j].unit_id != 65535)
			map2_usefulCount++;
		if ((*map1)[j].unit_id != 65535 && (*map2)[j].unit_id != 65535)
		{
			// printf("%s and %s\n", (*pl->names)[(*map1)[j].unit_id].c_str(),(*pl->names)[(*map2)[j].unit_id].c_str());
			// printf("%s and %s, %s\n", (*map1)[j].name.c_str(),(*map2)[j].name.c_str(), !((*map1)[j].forward ^ (*map2)[j].forward)?"forward":"backward");
			// printf("%s and %s\n", (*pl->names)[(*map1)[j].unit_id],(*pl->names)[(*map2)[j].unit_id]);
			// if(strcmp((*map1)[j].name,(*map2)[j].name)!=0){
			// 	failed_matches.push_back(j);
			// }
			int uid1 = (*map1)[j].unit_id;
			int uid2 = (*map2)[j].unit_id;
			int gid1 = 0, gid2 = 0;
			if (node_type_map1.find(uid1) != node_type_map1.end())
				gid1 = node_type_map1.at(uid1);
			if (node_type_map1.find(uid2) != node_type_map1.end())
				gid2 = node_type_map1.at(uid2);

			if (gid1 != 0 && gid2 != 0 && gid1 != gid2)
			{
				// 两个都不为0且不相等 → 跳过
				++skipped_count;
				continue;
			}
			success_counter_result++;
			// if(max_count ^ coverage[(*map1)[j].unit_id] !=0){
			coverage[(*map1)[j].unit_id]++;
			// }
			// if(max_count ^ coverage[(*map2)[j].unit_id] !=0){
			coverage[(*map2)[j].unit_id]++;
			if ((*map1)[j].unit_id != (*map2)[j].unit_id)
			{
				real_useful++;
			}
			// }
			// if(max_count ^ connections[(*map1)[j].unit_id][(*map2)[j].unit_id] !=0){
			if (!((*map1)[j].forward ^ (*map2)[j].forward))
			{
				connections_forward[(*map1)[j].unit_id][(*map2)[j].unit_id]++;
			}
			else
			{
				connections_backward[(*map1)[j].unit_id][(*map2)[j].unit_id]++;
			}
			// }

			std::string name1((*pl->names)[(*map1)[j].unit_id]);
			std::string name2((*pl->names)[(*map2)[j].unit_id]);

			// if(name1 != name2){
			// 	outFile << name1 << "," << (*map1)[j].pos << ","
            //     << name2 << "," << (*map2)[j].pos << "\n";
			// }
			
		}
	}

	printf("Useful Records: %ld\n", success_counter_result);
	printf("Not matched due to Wrong names: %ld\n", failed_matches.size());
	printf("real Useful Records: %ld\n", real_useful);
	printf("Skipped due to different group ids: %d\n", skipped_count);
	printf("Map1 count: %ld   Map2 count: %ld\n", map1_usefulCount, map2_usefulCount);
	if (out_fn)
	{
		FILE *fp = fopen(out_fn, "w");

		// for(int i = 0; i < pl->global_counter; i++){
		// 	uint32_t is_forward = forward[i] > backward[i] ? 1 : 0;
		// 	fprintf(fp, "%d\t%d\t%d\n", is_forward, backward[i], forward[i]);
		// }

		for (int i = 0; i < pl->global_counter; i++)
		{
			for (int j = 0; j < pl->global_counter; j++)
			{
				if (connections_forward[i][j] > 0 || connections_backward[i][j] > 0)
				{
					if (i < j)
					{
						connections_forward[j][i] += connections_forward[i][j];
						connections_backward[j][i] += connections_backward[i][j];
					}
					else if (i > j)
					{
						fprintf(fp, "%d\t%d\t%d\t%d\n", i, j, connections_forward[i][j], connections_backward[i][j]);
					}
				}
			}
		}
		fclose(fp);
	}

	// return connections;
}
void do_mapping(pldat_t *pl, bseq_file_t *hic_fn1, bseq_file_t *hic_fn2, char *out_fn)
{
	if (hic_fn1 == 0 || hic_fn2 == 0)
	{
		fprintf(stderr, "ERROR: Please give two hic files\n");
		exit(1);
	}

	ketopt_t o = KETOPT_INIT;
	int i, c, min_cnt = 2, mid_cnt = 5;
	tb_shared_t aux;

	memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = pl->opt->n_thread, aux.print_diff = 0;
	aux.ratio_thres = 0.33;
	aux.k = pl->h->k;
	aux.ch = pl->h;
	aux.h_pos = pl->h_pos; // 新增
	aux.record_num = pl->global_counter;
	cout << "Start mapping, total unitig number: " << aux.record_num << endl;
	aux.fp = hic_fn1;
	std::vector<mapping_res_t> *map1 = new std::vector<mapping_res_t>();
	aux.mappings = map1;
	aux.buf = (tb_buf_t *)calloc(aux.n_threads, sizeof(tb_buf_t));
	cout << "Start hic1 mapping " << endl;
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	for (i = 0; i < aux.n_threads; ++i)
	{
		free(aux.buf[i].s);
	}
	free(aux.buf);

	aux.fp = hic_fn2;
	std::vector<mapping_res_t> *map2 = new std::vector<mapping_res_t>();
	aux.mappings = map2;
	aux.buf = (tb_buf_t *)calloc(aux.n_threads, sizeof(tb_buf_t));
	cout << "Start hic2 mapping " << endl;
	kt_pipeline(2, tb_pipeline, &aux, 2);
	bseq_close(aux.fp);
	for (i = 0; i < aux.n_threads; ++i)
	{
		free(aux.buf[i].s);
	}
	free(aux.buf);

	printf("Map1 size: %ld;\n", map1->size());
	printf("Map2 size: %ld;\n", map2->size());

	uint64_t success_counter_result = 0;
	uint32_t **connections_forward;
	CALLOC(connections_forward, pl->global_counter);
	for (int i = 0; i < pl->global_counter; i++)
	{
		CALLOC(connections_forward[i], pl->global_counter);
		memset(connections_forward[i], 0, sizeof(connections_forward[i]));
	}
	uint32_t **connections_backward;
	CALLOC(connections_backward, pl->global_counter);
	for (int i = 0; i < pl->global_counter; i++)
	{
		CALLOC(connections_backward[i], pl->global_counter);
		memset(connections_backward[i], 0, sizeof(connections_backward[i]));
	}
	uint32_t coverage[pl->global_counter];
	memset(coverage, 0, sizeof(uint32_t) * pl->global_counter);
	std::vector<uint64_t> failed_matches;
	const uint16_t max_count = -1;
	for (uint64_t j = 0; j < std::min(map1->size(), map2->size()); j++)
	{
		// if((*map1)[j].unit_id != 65535){
		// 	if((*map1)[j].forward){
		// 		forward[(*map1)[j].unit_id]++;
		// 	}else{
		// 		backward[(*map1)[j].unit_id]++;
		// 	}
		// }
		// if((*map2)[j].unit_id != 65535){
		// 	if((*map2)[j].forward){
		// 		forward[(*map2)[j].unit_id]++;
		// 	}else{
		// 		backward[(*map2)[j].unit_id]++;
		// 	}
		// }
		if ((*map1)[j].unit_id != 65535 && (*map2)[j].unit_id != 65535)
		{
			// printf("%s and %s\n", (*pl->names)[(*map1)[j].unit_id].c_str(),(*pl->names)[(*map2)[j].unit_id].c_str());
			// printf("%s and %s, %s\n", (*map1)[j].name.c_str(),(*map2)[j].name.c_str(), !((*map1)[j].forward ^ (*map2)[j].forward)?"forward":"backward");
			// printf("%s and %s\n", (*pl->names)[(*map1)[j].unit_id],(*pl->names)[(*map2)[j].unit_id]);
			// if(strcmp((*map1)[j].name,(*map2)[j].name)!=0){
			// 	failed_matches.push_back(j);
			// }

			success_counter_result++;
			// if(max_count ^ coverage[(*map1)[j].unit_id] !=0){
			coverage[(*map1)[j].unit_id]++;
			// }
			// if(max_count ^ coverage[(*map2)[j].unit_id] !=0){
			coverage[(*map2)[j].unit_id]++;
			// }
			// if(max_count ^ connections[(*map1)[j].unit_id][(*map2)[j].unit_id] !=0){
			if (!((*map1)[j].forward ^ (*map2)[j].forward))
			{
				connections_forward[(*map1)[j].unit_id][(*map2)[j].unit_id]++;
			}
			else
			{
				connections_backward[(*map1)[j].unit_id][(*map2)[j].unit_id]++;
			}
			// }
		}
	}

	printf("Useful Records: %ld\n", success_counter_result);
	printf("Not matched due to Wrong names: %ld\n", failed_matches.size());

	if (out_fn)
	{
		FILE *fp = fopen(out_fn, "w");

		// for(int i = 0; i < pl->global_counter; i++){
		// 	uint32_t is_forward = forward[i] > backward[i] ? 1 : 0;
		// 	fprintf(fp, "%d\t%d\t%d\n", is_forward, backward[i], forward[i]);
		// }

		for (int i = 0; i < pl->global_counter; i++)
		{
			for (int j = 0; j < pl->global_counter; j++)
			{
				if (connections_forward[i][j] > 0 || connections_backward[i][j] > 0)
				{
					if (i < j)
					{
						connections_forward[j][i] += connections_forward[i][j];
						connections_backward[j][i] += connections_backward[i][j];
					}
					else if (i > j)
					{
						fprintf(fp, "%d\t%d\t%d\t%d\n", i, j, connections_forward[i][j], connections_backward[i][j]);
					}
				}
			}
		}
		fclose(fp);
	}

	// return connections;
}

void load_node_type_csv(const char *csv_filename, std::unordered_map<int, int> &node_type_map)
{
	std::ifstream fin(csv_filename);
	if (!fin.is_open())
	{
		fprintf(stderr, "ERROR: cannot open CSV file %s\n", csv_filename);
		exit(1);
	}

	std::string line;
	bool first_line = true;
	while (std::getline(fin, line))
	{
		if (line.empty())
			continue;

		// 如果CSV第一行是表头，跳过
		if (first_line && (line.find("node") != std::string::npos || line.find("Node") != std::string::npos))
		{
			first_line = false;
			continue;
		}
		first_line = false;

		std::stringstream ss(line);
		std::string id_str, type_str;

		if (std::getline(ss, id_str, ',') && std::getline(ss, type_str))
		{
			try
			{
				int node_id = std::stoi(id_str);
				int node_type = std::stoi(type_str);
				node_type_map[node_id] = node_type;
			}
			catch (const std::exception &e)
			{
				std::cerr << "WARNING: skip invalid line in CSV: " << line << std::endl;
			}
		}
	}
	fin.close();
}

int main_hic_map(int argc, char *argv[])
{
	pldat_t *h;
	int c;
	char *fn_out = 0;
	yak_copt_t opt;
	ketopt_t o = KETOPT_INIT;
	yak_copt_init(&opt);
	opt.pre = YAK_COUNTER_BITS;
	opt.n_thread = 32;
	while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:o:", 0)) >= 0)
	{
		if (c == 'k')
			opt.k = atoi(o.arg);
		else if (c == 'p')
			opt.pre = atoi(o.arg);
		else if (c == 'K')
			opt.chunk_size = atoi(o.arg);
		else if (c == 't')
			opt.n_thread = atoi(o.arg);
		else if (c == 'b')
			opt.bf_shift = atoi(o.arg);
		else if (c == 'o')
			fn_out = o.arg;
	}
	if (argc - o.ind < 1)
	{
		fprintf(stderr, "Usage: pstools hic_mapping [options] <graph.fa> <hic.R1.fastq.gz> <hic.R2.fastq.gz>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", opt.pre);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", opt.bf_shift);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "  -o FILE    save mapping relationship to FILE []\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		fprintf(stderr, "Test\n");
		return 1;
	}
	if (opt.pre < YAK_COUNTER_BITS)
	{
		fprintf(stderr, "ERROR: -p should be at least %d\n", YAK_COUNTER_BITS);
		return 1;
	}
	if (opt.k >= 64)
	{
		fprintf(stderr, "ERROR: -k must be smaller than 64\n");
		return 1;
	}
	else if (opt.k >= 32)
	{
		fprintf(stderr, "WARNING: counts are inexact if -k is greater than 31\n");
	}

	if (argv[o.ind + 1] == 0 || argv[o.ind + 2] == 0)
	{
		fprintf(stderr, "ERROR: Please give two hic files\n");
		exit(1);
	}

	bseq_file_t *hic_fn1 = bseq_open(argv[o.ind + 1]);
	bseq_file_t *hic_fn2 = bseq_open(argv[o.ind + 2]);
	if (hic_fn1 == 0 || hic_fn2 == 0)
	{
		fprintf(stderr, "ERROR: fail to open hic files\n");
		exit(1);
	}

	h = yak_count_multi_new(argv[o.ind], &opt, 0);

	do_mapping(h, hic_fn1, hic_fn2, fn_out);

	free(h->h);
	free(h);
	return 0;
}

//......................................polishing.................................................

static void hifi_worker(void *_data, long k, int tid)
{
	tb_step_hifi *t = (tb_step_hifi *)_data;
	tb_shared_hifi_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	tb_buf_t *b = &aux->buf[tid];
	uint64_t x[4], mask;
	int i, l, shift;
	bool is_forward;
	std::map<uint16_t, int> counts;
	std::map<uint16_t, bool> forward;

	std::map<uint16_t, int> counts_for_kmer_h1;
	std::map<uint16_t, bool> forward_for_kmer_h1;
	std::map<uint16_t, int> counts_for_kmer_h2;
	std::map<uint16_t, bool> forward_for_kmer_h2;
	// Porec
	std::vector<Hamper_kmer> hamper_kmer;
	std::map<int, std::vector<Tread2shapmers>> read2hapmers;
	std::map<int, std::vector<Tread2shapmers>> read2hapmers_hap1;
	std::map<int, std::vector<Tread2shapmers>> read2hapmers_hap2;
	std::unordered_map<uint64_t, int> hap_length_map;
	// 。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。
	bool hap_flag = false;
	if (aux->ch1->k < 32)
	{
		mask = (1ULL << 2 * aux->ch1->k) - 1;
		shift = 2 * (aux->ch1->k - 1);
	}
	else
	{
		mask = (1ULL << aux->ch1->k) - 1;
		shift = aux->ch1->k - 1;
	}
	if (s->l_seq > b->max)
	{
		b->max = s->l_seq;
		kroundup32(b->max);
		b->s = (uint32_t *)realloc(b->s, b->max * sizeof(uint32_t));
	}
	memset(b->s, 0, s->l_seq * sizeof(uint32_t));
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->l_seq; ++i)
	{
		int c = seq_nt4_table[(uint8_t)s->seq[i]];
		char qual = s->qual[i];
		int res_h1, res_h2;
		int res_h1_kmer, res_h2_kmer;
		int res_h1_pos, res_h2_pos;
		if (c < 4)
		{
			if (aux->ch1->k < 32)
			{
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			}
			else
			{
				x[0] = (x[0] << 1 | (c & 1)) & mask;
				x[1] = (x[1] << 1 | (c >> 1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c & 1)) << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c >> 1)) << shift;
			}
			if (++l >= aux->ch1->k)
			{

				uint64_t y, y1; // y1 for kmer,y for hamper

				if (aux->ch1->k < 32)
				{
					if (qual == 'H')
					{
						is_forward = x[0] < x[1];
						y = yak_hash64(is_forward ? x[0] : x[1], mask);
					}
					else
					{
						y1 = yak_hash64(is_forward ? x[0] : x[1], mask);
					}
				}
				else
				{
					y = yak_hash_long(x);
				}

				res_h1 = yak_ch_get_k(aux->ch1, y);
				res_h2 = yak_ch_get_k(aux->ch2, y);
				res_h1_kmer = yak_ch_get_k(aux->ch1, y1);
				res_h2_kmer = yak_ch_get_k(aux->ch2, y1);

				if (res_h1_kmer != -1)
				{
					counts_for_kmer_h1[res_h1_kmer & YAK_KEY_MASK]++;
					bool is_main_forward = !(((res_h1_kmer & YAK_FORWARD_MASK) != 0) ^ is_forward);
					forward_for_kmer_h1[res_h1_kmer & YAK_KEY_MASK] = is_main_forward;
					// cout<<"res_h1_kmer: "<<res_h1_kmer<<" "<<is_main_forward<<" "<<res_h1_kmer<<" "<<(res_h1_kmer & YAK_KEY_MASK)<<" "<<counts_for_kmer_h1[res_h1_kmer & YAK_KEY_MASK]<<" "<<res_h1_kmer<<" "<<(res_h1_kmer & YAK_KEY_MASK)<<endl;
				}
				if (res_h2_kmer != -1)
				{
					counts_for_kmer_h2[res_h2_kmer & YAK_KEY_MASK]++;
					bool is_main_forward = !(((res_h2_kmer & YAK_FORWARD_MASK) != 0) ^ is_forward);
					forward_for_kmer_h2[res_h2_kmer & YAK_KEY_MASK] = is_main_forward;
				}

				if (res_h1 != -1 || res_h2 != -1)
				{
					if (res_h1 != -1)
					{

						uint32_t pos = -1;
						uint16_t read_id = yak_ch_get_pos(aux->ch1, aux->h_pos1, y, &pos);
						if (pos != -1)
						{
							int target_id = res_h1 & YAK_KEY_MASK;
							Tread2shapmers hap;
							// hap.hapmer = decode_kmer(x[0], aux->ch1->k);
							hap.hapmer = x[0];

							hap.tposi = pos;
							hap.strand = !(((res_h1 & YAK_FORWARD_MASK) != 0) ^ is_forward); // 相同为0，不同为1

							hap.qposi = l - aux->ch1->k;

							read2hapmers_hap1[target_id].push_back(hap);
						}
					}

					if (res_h2 != -1)
					{
						uint32_t pos = -1;
						uint16_t read_id = yak_ch_get_pos(aux->ch2, aux->h_pos2, y, &pos);
						if (pos != -1)
						{
							int target_id = res_h2 & YAK_KEY_MASK;
							Tread2shapmers hap;
							// hap.hapmer = decode_kmer(x[0], aux->ch2->k);
							hap.hapmer = x[0];

							hap.tposi = pos;
							hap.strand = !(((res_h2 & YAK_FORWARD_MASK) != 0) ^ is_forward); // 相同为0，不同为1

							hap.qposi = l - aux->ch2->k;

							read2hapmers_hap2[target_id].push_back(hap);
						}
					}
				}
			}
		}
		else
			l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}

	// TODO:
	if (read2hapmers_hap1.size() == 0 && read2hapmers_hap2.size() == 0)
	{
		return;
	}

	// 替换原有的判断块

	// 找出 map 中 value.size() 的最大值（如果 map 为空，返回 0）
	auto max_vec_size = [](const auto &m) -> size_t
	{
		size_t mx = 0;
		for (const auto &kv : m)
		{
			mx = std::max(mx, kv.second.size());
		}
		return mx;
	};

	size_t max1 = max_vec_size(read2hapmers_hap1);
	size_t max2 = max_vec_size(read2hapmers_hap2);

	if (max1 > max2)
	{
		read2hapmers = read2hapmers_hap1;
		hap_length_map = g_read_length_map_hap1;
		hap_flag = true;
	}
	else if (max2 > max1)
	{
		read2hapmers = read2hapmers_hap2;
		hap_length_map = g_read_length_map_hap2;
		hap_flag = false;
	}
	else
	{
		// tie-break: 可以按你希望的策略决定，这里示例保留原逻辑改成 hap1 优先
		read2hapmers = read2hapmers_hap1;
		hap_length_map = g_read_length_map_hap1;
		hap_flag = true;
	}

	// get max kmer number and keys
	uint64_t max_key_h1 = 0;
	int max_value_h1 = -1;

	for (const auto &p : counts_for_kmer_h1)
	{
		if (p.second > max_value_h1)
		{
			max_value_h1 = p.second;
			max_key_h1 = p.first;
		}
	}
	// cout<<"max_key_h1: "<<max_key_h1<<" "<<max_value_h1<<endl;

	uint64_t max_key_h2 = 0;
	int max_value_h2 = -1;

	for (const auto &p : counts_for_kmer_h2)
	{
		if (p.second > max_value_h2)
		{
			max_value_h2 = p.second;
			max_key_h2 = p.first;
		}
	}

	for (const auto &pair : read2hapmers)
	{
		int target_id = pair.first;
		const std::vector<Tread2shapmers> &shapmers = pair.second;
		int num_shapmers = shapmers.size();

		// std::cout << "Target ID: " << target_id << "\n";
		// for (const auto &h : shapmers)
		// {
		// 	std::cout << "  Hapmer: " << h.hapmer
		// 			  << ", QPos: " << h.qposi
		// 			  << ", TPos: " << h.tposi
		// 			  << ", Strand: " << h.strand << "\n";
		// }

		if (num_shapmers < min_num_shapmers)
		{
			continue;
		}
		// 计算相同 strand 和不同 strand 的数量
		int num_same_strand = 0;
		int num_diff_strand = 0;

		// 遍历 shapmers，统计相同和不同 strand 的数量
		for (const Tread2shapmers &shapmer : shapmers)
		{
			if (shapmer.strand == 1)
			{
				num_same_strand++;
			}
			else
			{
				num_diff_strand++;
			}
		}
		// 确定是同向比对，还是反向互补的比对上
		int strand = (num_same_strand > num_diff_strand) ? 1 : 0;

		std::vector<Tread2shapmers> strand_correct_shapmers;
		// 遍历 shapmers，只挑出正确比对方向的hamper，过滤掉错误噪音
		for (const Tread2shapmers &shapmer : shapmers)
		{
			if (shapmer.strand == strand)
			{
				strand_correct_shapmers.push_back(shapmer);
			}
		}

		// std::cout << "	Target ID in strand_correct_shapmers: " << target_id << "\n";
		// for (const auto &h : strand_correct_shapmers)
		// {
		// 	std::cout << "  	Hapmer: " << h.hapmer
		// 			  << "	, QPos: " << h.qposi
		// 			  << "	, TPos: " << h.tposi
		// 			  << "	, Strand: " << h.strand << "\n";
		// }

		// if (strand == 0)
		// {

		// }

		// std::sort(strand_correct_shapmers.begin(), strand_correct_shapmers.end(),
		// 		  [](const Tread2shapmers &a, const Tread2shapmers &b)
		// 		  {
		// 			  return a.qposi < b.qposi;
		// 		  });
		// 统计符合 strand 的 hapmers 数量
		int num_strand_correct_shapmers = strand_correct_shapmers.size();
		if (num_strand_correct_shapmers < min_num_shapmers)
		{
			continue;
		}
		std::vector<Tread2shapmers> posi_correct_shapmers;
		int n = 3;		 // 使用 n 个下一个 hapmers 来判断当前 hapmer
		int shift = 200; // 允许的最大位置偏移
		if (num_strand_correct_shapmers >= 2 * n)
		{
			if (strand == 1) // 相同的 strand
			{
				for (int i = 0; i <= num_strand_correct_shapmers - n - 1; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i + j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) <= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
						{
							score += -1;
						}
						else
						{
							score += 1;
						}
					}

					if (score > 0) // 认为这个 hapmer 是可信的，否则移除
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}

				// 处理最后几个 hapmers（这些 hapmers 没有 n 个后续 hapmers 来判断）
				for (int i = num_strand_correct_shapmers - n; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i - j];
						int qposi_p = hapmer_j.qposi;
						int tposi_p = hapmer_j.tposi;

						if ((qposi - qposi_p) <= 0 || (tposi - tposi_p) <= 0 || std::abs(std::abs(qposi - qposi_p) - std::abs(tposi - tposi_p)) > shift)
						{
							score += -1;
						}
						else
						{
							score += 1;
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
			//
			else // 不相同的 strand
			{
				for (int i = 0; i <= num_strand_correct_shapmers - n - 1; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i + j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) >= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
						{
							score += -1;
							// cout<<"Now score -1 : "<<score<<endl;
						}
						else
						{
							score += 1;
							// cout<<"Now score +1 : "<<score<<endl;
						}
					}

					if (score > 0) // 认为这个 hapmer 是可信的，否则移除
					{
						// cout<<"Accepted hapmer: "<<strand_correct_shapmers[i].hapmer<<", QPos: "<<strand_correct_shapmers[i].qposi<<", TPos: "<<strand_correct_shapmers[i].tposi<<", Strand: "<<strand_correct_shapmers[i].strand<<endl;
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}

				// 处理最后几个 hapmers（这些 hapmers 没有 n 个后续 hapmers 来判断）
				for (int i = num_strand_correct_shapmers - n; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i - j];
						int qposi_p = hapmer_j.qposi;
						int tposi_p = hapmer_j.tposi;

						if ((qposi - qposi_p) <= 0 || (tposi - tposi_p) >= 0 || std::abs(std::abs(qposi - qposi_p) - std::abs(tposi - tposi_p)) > shift)
						{
							score += -1;
						}
						else
						{
							score += 1;
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
		}
		else if (num_strand_correct_shapmers >= 2 && num_strand_correct_shapmers < 2 * n)
		{
			if (strand == 1) // 相同的 strand
			{
				for (int i = 0; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 0; j < num_strand_correct_shapmers; ++j)
					{
						if (i == j)
						{
							continue;
						}

						const auto &hapmer_j = strand_correct_shapmers[j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if (j > i)
						{
							if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) <= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
						else
						{
							if ((qposi_n - qposi) >= 0 || (tposi_n - tposi) >= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
			else // 不同的 strand
			{
				for (int i = 0; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 0; j < num_strand_correct_shapmers; ++j)
					{
						if (i == j)
						{
							continue;
						}

						const auto &hapmer_j = strand_correct_shapmers[j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if (j > i)
						{
							if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) >= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
						else
						{
							if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) <= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
		}
		else
		{
			continue;
		}
		if (posi_correct_shapmers.size() < min_num_shapmers)
		{
			continue;
		}

		int num_nonovlp_shapmers = 1; // 初始值为 1，考虑第一个 hapmer

		std::unordered_map<uint64_t, int> count;
		for (const auto &h : posi_correct_shapmers)
		{
			count[h.hapmer]++;
		}

		std::vector<Tread2shapmers> filtered;
		for (const auto &h : posi_correct_shapmers)
		{
			if (count[h.hapmer] == 1)
			{ // 只保留出现 1 次的
				filtered.push_back(h);
			}
		}
		posi_correct_shapmers.swap(filtered);
		if (filtered.size() < min_num_shapmers)
		{
			continue;
		}
		cout << "Quary name:" << s->name << endl;
		for (size_t i = 0; i < filtered.size(); ++i)
		{
			std::cout << "After filtering, Hapmer: " << filtered[i].hapmer
					  << "	, QPos: " << filtered[i].qposi
					  << "	, TPos: " << filtered[i].tposi
					  << "	, Strand: " << filtered[i].strand << "\n";
		}

		//...............................................................................

		for (size_t i = 0; i < filtered.size(); ++i)
		{
			const Tread2shapmers &hapmer1 = filtered[i];
			const Tread2shapmers &hapmer2 = filtered[i + 1];

			int qposi1 = hapmer1.qposi;
			int tposi1 = hapmer1.tposi;
			int qposi2 = hapmer2.qposi;
			int tposi2 = hapmer2.tposi;
			// cout << "  Filtered Hapmer: " << hapmer1.hapmer
			// 	 << ", QPos: " << hapmer1.qposi
			// 	 << ", TPos: " << hapmer1.tposi
			// 	 << ", Strand: " << hapmer1.strand << "\n";

			// 检查相邻 hapmers 是否满足距离要求
			if (std::abs(qposi2 - qposi1) >= aux->ch1->k && std::abs(tposi2 - tposi1) >= aux->ch1->k)
			{
				num_nonovlp_shapmers++;
			}
		}

		if (posi_correct_shapmers.size() > 60)
		{
			int num = posi_correct_shapmers.size() / 60;
			num_nonovlp_shapmers += num;
		}

		const Tread2shapmers &first_hapmer = filtered.front();
		const Tread2shapmers &last_hapmer = filtered.back();
		int qposi1 = first_hapmer.qposi;
		int tposi1 = first_hapmer.tposi;
		int qposi2 = last_hapmer.qposi;
		int tposi2 = last_hapmer.tposi;
		strand = first_hapmer.strand;
		int qstart, qend, tstart, tend;

		if (strand == 1)
		{ // same strand
			qstart = qposi1;
			qend = qposi2;
			tstart = tposi1;
			tend = tposi2;
		}
		else
		{ // different strand
			qstart = qposi1;
			qend = qposi2;
			tstart = tposi2;
			tend = tposi1;
		}
		// if (min(qend-qstart, tend-tstart)/max(qend-qstart, tend-tstart)) < 0.5:
		//     continue
		// // // 检查重叠比例
		double overlap = static_cast<double>(std::min(qend - qstart, tend - tstart)) / std::max(qend - qstart, tend - tstart);
		// cout << "	overlap ratio: " << overlap << " " << static_cast<double>(std::min(qend - qstart, tend - tstart)) << " " << std::max(qend - qstart, tend - tstart) << endl;
		if (num_nonovlp_shapmers < 2)
		{
			continue;
		}

		if (static_cast<double>(std::min(qend - qstart, tend - tstart)) / std::max(qend - qstart, tend - tstart) < 0.5)
		{
			continue;
		}
		int qspan = std::abs(qend - qstart);
		int tspan = std::abs(tend - tstart);
		int min_len = hap_length_map[target_id] < s->l_seq ? hap_length_map[target_id] : s->l_seq;

		// qspan < min_len / 6 or tspan < min_len / 6
		if (qspan < min_len / 10 || tspan < min_len / 10)
		{
			continue;
		}

		// cout << "Num non-overlapping shapmers: " << num_nonovlp_shapmers << "\n";
		counts[target_id] = num_nonovlp_shapmers;
		forward[target_id] = strand;
		if (hap_flag == true)
		{
			// cout<<"Selected hap1 target_id: "<<target_id <<",max_key_h1:"<<max_key_h1<<endl;
			if (max_key_h1 != target_id)
			{
				continue;
			}
			std::ofstream outFile;
			outFile.open("output_hap1.csv", std::ios::app);
			outFile << (*aux->contigs_names_hap1)[target_id] << ","
					<< hap_length_map[target_id] << ","
					<< tstart << ","
					<< tend << ","
					<< strand << ","
					<< s->name << ","
					<< s->l_seq << ","
					<< qstart << ","
					<< qend << ","
					<< num_nonovlp_shapmers << ","
					<< filtered.size() << ","
					<< posi_correct_shapmers.size() << "\n";
			// 关闭文件
			outFile.close();
		}
		else
		{

			if (max_key_h2 != target_id)
			{
				continue;
			}
			std::ofstream outFile;
			outFile.open("output_hap2.csv", std::ios::app);
			outFile << (*aux->contigs_names_hap2)[target_id] << ","
					<< hap_length_map[target_id] << ","
					<< tstart << ","
					<< tend << ","
					<< strand << ","
					<< s->name << ","
					<< s->l_seq << ","
					<< qstart << ","
					<< qend << ","
					<< num_nonovlp_shapmers << ","
					<< filtered.size() << ","
					<< posi_correct_shapmers.size() << "\n";
			// 关闭文件
			outFile.close();
		}

		if (hap_flag == true)
		{
			std::ofstream pafFile;
			pafFile.open("output_hap1.paf", std::ios::app);

			// 提取字段
			const std::string &target_name = (*aux->contigs_names_hap1)[target_id];
			int tlen = hap_length_map[target_id];
			int qlen = s->l_seq;
			int qspan = std::abs(qend - qstart);
			int tspan = std::abs(tend - tstart);
			int aln_block_len = std::max(qspan, tspan);
			double identity = 0.85;
			int n_match = static_cast<int>(aln_block_len * identity);

			// MAPQ 计算
			int mapq = 0;
			if (num_nonovlp_shapmers < 5)
				mapq = 0;
			else if (num_nonovlp_shapmers > 20)
				mapq = 60;
			else
			{
				mapq = static_cast<int>(1 + (num_nonovlp_shapmers - 5) * (58.0 / 25.0));
				if (mapq > 59)
					mapq = 59;
			}

			// strand 转换为 '+' 或 '-'
			char strand_symbol = (strand == 1 ? '+' : '-');

			// 写入 PAF 行
			pafFile << s->name << "\t"
					<< qlen << "\t"
					<< qstart << "\t"
					<< qend << "\t"
					<< strand_symbol << "\t"
					<< target_name << "\t"
					<< tlen << "\t"
					<< tstart << "\t"
					<< tend << "\t"
					<< n_match << "\t"
					<< aln_block_len << "\t"
					<< mapq << "\n";

			pafFile.close();
		}
		else
		{
			std::ofstream pafFile;
			pafFile.open("output_hap2.paf", std::ios::app);

			const std::string &target_name = (*aux->contigs_names_hap2)[target_id];
			int tlen = hap_length_map[target_id];
			int qlen = s->l_seq;
			int qspan = std::abs(qend - qstart);
			int tspan = std::abs(tend - tstart);
			int aln_block_len = std::max(qspan, tspan);
			double identity = 0.85;
			int n_match = static_cast<int>(aln_block_len * identity);

			// MAPQ 计算
			int mapq = 0;
			if (num_nonovlp_shapmers < 5)
				mapq = 0;
			else if (num_nonovlp_shapmers > 20)
				mapq = 60;
			else
			{
				mapq = static_cast<int>(1 + (num_nonovlp_shapmers - 5) * (58.0 / 25.0));
				if (mapq > 59)
					mapq = 59;
			}

			char strand_symbol = (strand == 1 ? '+' : '-');

			pafFile << s->name << "\t"
					<< qlen << "\t"
					<< qstart << "\t"
					<< qend << "\t"
					<< strand_symbol << "\t"
					<< target_name << "\t"
					<< tlen << "\t"
					<< tstart << "\t"
					<< tend << "\t"
					<< n_match << "\t"
					<< aln_block_len << "\t"
					<< mapq << "\n";

			pafFile.close();
		}
	}
}

static void porec_worker(void *_data, long k, int tid)
{
	tb_step_t *t = (tb_step_t *)_data;
	tb_shared_t *aux = t->aux;
	bseq1_t *s = &t->seq[k];
	tb_buf_t *b = &aux->buf[tid];
	uint64_t x[4], mask;
	int i, l, shift;
	bool is_forward;
	std::map<uint16_t, int> counts;
	std::map<uint16_t, bool> forward;

	// Porec
	std::vector<Hamper_kmer> hamper_kmer;
	std::map<int, std::vector<Tread2shapmers>> read2hapmers;
	// 。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。

	if (aux->ch->k < 32)
	{
		mask = (1ULL << 2 * aux->ch->k) - 1;
		shift = 2 * (aux->ch->k - 1);
	}
	else
	{
		mask = (1ULL << aux->ch->k) - 1;
		shift = aux->ch->k - 1;
	}
	if (s->l_seq > b->max)
	{
		b->max = s->l_seq;
		kroundup32(b->max);
		b->s = (uint32_t *)realloc(b->s, b->max * sizeof(uint32_t));
	}
	memset(b->s, 0, s->l_seq * sizeof(uint32_t));
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->l_seq; ++i)
	{
		int c = seq_nt4_table[(uint8_t)s->seq[i]];
		int res;
		if (c < 4)
		{
			if (aux->ch->k < 32)
			{
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
			}
			else
			{
				x[0] = (x[0] << 1 | (c & 1)) & mask;
				x[1] = (x[1] << 1 | (c >> 1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c & 1)) << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c >> 1)) << shift;
			}
			if (++l >= aux->ch->k)
			{

				uint64_t y;

				if (aux->ch->k < 32)
				{
					is_forward = x[0] < x[1];
					y = yak_hash64(is_forward ? x[0] : x[1], mask);
				}
				else
				{
					y = yak_hash_long(x);
				}
				res = yak_ch_get_k(aux->ch, y);
				int res1 = yak_ch_get_k_for_repeat(aux->ch, y);
				// cout << "kmer: " << decode_kmer(is_forward ? x[0] : x[1], aux->ch->k) << "is_forward:" << is_forward << ",x[0] kmer: " << decode_kmer(x[0], aux->ch->k)
				// 	 << ", x[1] kmer: " << decode_kmer(x[1], aux->ch->k) << ", y: " << y << ", res: " << res << ", res1: " << res1 << "\n";
				if (res1 != -1 && res == -1)
				{

					//  检查是否重复
					auto it = repeated_kmers_map.find(y);
					if (it != repeated_kmers_map.end())
					{
						// cout << "y=" << y << " is a repeated kmer, repeat times: " << it->second.size() << "\n";
						for (auto &info : it->second)
						{

							// int target_id = res1 & YAK_KEY_MASK;
							Tread2shapmers hap;
							// hap.hapmer = decode_kmer(x[0], aux->ch->k);
							hap.hapmer = x[0];

							hap.tposi = info.pos;
							hap.strand = !(info.forward ^ is_forward); // 相同为0，不同为1

							hap.qposi = l - aux->ch->k;

							read2hapmers[info.read_name].push_back(hap);
							// if (info.read_name == 0)
							// {
							// 	cout << "Other hampers in vector" << ",hapmer: " << hap.hapmer << ", qposi: " << hap.qposi << ", tposi: " << hap.tposi << ", strand: " << hap.strand << "\n";
							// }

							// cout << "info.read_name:" << info.read_name << ",hamper: " << hap.hapmer << ", qposi: " << hap.qposi << ", tposi: " << hap.tposi << ", strand: " << hap.strand << "\n";
						}
					}
					else
					{
						cout << "y=" << y << " is a new repeated kmer\n";
					}
					// TODO: 对于第一次出现的 hamper，选择记录
					uint32_t pos = -1;
					uint64_t read_id = yak_ch_get_pos_for_repeat(aux->ch, aux->h_pos, y, &pos);
					if (pos != -1)
					{
						int target_id = res1 & YAK_KEY_MASK;
						Tread2shapmers hap;
						// hap.hapmer = decode_kmer(x[0], aux->ch->k);
						hap.hapmer = x[0];
						hap.tposi = pos;
						hap.strand = !(((res1 & YAK_FORWARD_MASK) != 0) ^ is_forward); // 相同为0，不同为1

						hap.qposi = l - aux->ch->k;

						read2hapmers[target_id].push_back(hap);
						// if (target_id == 0)
						// {
						// 	cout << "First  hamper in vector" << ",hapmer: " << hap.hapmer << ", qposi: " << hap.qposi << ", tposi: " << hap.tposi << ", strand: " << hap.strand << "\n";
						// }
						// cout << "target_id:" << target_id << ",hapmer: " << hap.hapmer << ", qposi: " << hap.qposi << ", tposi: " << hap.tposi << ", strand: " << hap.strand << "\n";
					}
				}
				if (res != -1)
				{
					// cout << "is_forward:" << is_forward << ",x[0] kmer: " << decode_kmer(x[0], aux->ch->k) << ", x[1] kmer: " << decode_kmer(x[1], aux->ch->k) << "\n";
					uint32_t pos = -1;
					uint16_t read_id = yak_ch_get_pos(aux->ch, aux->h_pos, y, &pos);
					if (pos != -1)
					{
						int target_id = res & YAK_KEY_MASK;
						Tread2shapmers hap;
						//hap.hapmer = decode_kmer(x[0], aux->ch->k);
						hap.hapmer = x[0];
						hap.tposi = pos;
						hap.strand = !(((res & YAK_FORWARD_MASK) != 0) ^ is_forward); // 相同为0，不同为1

						hap.qposi = l - aux->ch->k;

						read2hapmers[target_id].push_back(hap);
						// if (target_id == 0)
						// {
						// 	cout << "only hamper in hash table" << ",hapmer: " << hap.hapmer << ", qposi: " << hap.qposi << ", tposi: " << hap.tposi << ", strand: " << hap.strand << "\n";
						// }
						// cout << "y=" << y << ", read_id=" << read_id << ", pos=" << pos << ", target_id=" << target_id << ", strand=" << hap.strand << "\n";
						// cout << "first target_id: " << target_id << ",hapmer: " << hap.hapmer << ", qposi: " << hap.qposi << ", tposi: " << hap.tposi << ", strand: " << hap.strand << "\n";
					}
				}
			}
		}
		else
			l = 0, x[0] = x[1] = x[2] = x[3] = 0;
	}

	// TODO:

	for (const auto &pair : read2hapmers)
	{
		int target_id = pair.first;
		const std::vector<Tread2shapmers> &shapmers = pair.second;
		int num_shapmers = shapmers.size();

		// std::cout << "Target ID: " << target_id << "\n";
		// for (const auto &h : shapmers)
		// {
		// 	std::cout << "  Hapmer: " << h.hapmer
		// 			  << ", QPos: " << h.qposi
		// 			  << ", TPos: " << h.tposi
		// 			  << ", Strand: " << h.strand << "\n";
		// }

		if (num_shapmers < min_num_shapmers)
		{
			continue;
		}
		// 计算相同 strand 和不同 strand 的数量
		int num_same_strand = 0;
		int num_diff_strand = 0;

		// 遍历 shapmers，统计相同和不同 strand 的数量
		for (const Tread2shapmers &shapmer : shapmers)
		{
			if (shapmer.strand == 1)
			{
				num_same_strand++;
			}
			else
			{
				num_diff_strand++;
			}
		}
		// 确定是同向比对，还是反向互补的比对上
		int strand = (num_same_strand > num_diff_strand) ? 1 : 0;

		std::vector<Tread2shapmers> strand_correct_shapmers;
		// 遍历 shapmers，只挑出正确比对方向的hamper，过滤掉错误噪音
		for (const Tread2shapmers &shapmer : shapmers)
		{
			if (shapmer.strand == strand)
			{
				strand_correct_shapmers.push_back(shapmer);
			}
		}

		// std::cout << "	Target ID in strand_correct_shapmers: " << target_id << "\n";
		// for (const auto &h : strand_correct_shapmers)
		// {
		// 	std::cout << "  	Hapmer: " << h.hapmer
		// 			  << "	, QPos: " << h.qposi
		// 			  << "	, TPos: " << h.tposi
		// 			  << "	, Strand: " << h.strand << "\n";
		// }

		// if (strand == 0)
		// {

		// }

		// std::sort(strand_correct_shapmers.begin(), strand_correct_shapmers.end(),
		// 		  [](const Tread2shapmers &a, const Tread2shapmers &b)
		// 		  {
		// 			  return a.qposi < b.qposi;
		// 		  });
		// 统计符合 strand 的 hapmers 数量
		int num_strand_correct_shapmers = strand_correct_shapmers.size();
		if (num_strand_correct_shapmers < min_num_shapmers)
		{
			continue;
		}
		std::vector<Tread2shapmers> posi_correct_shapmers;
		int n = 3;		 // 使用 n 个下一个 hapmers 来判断当前 hapmer
		int shift = 200; // 允许的最大位置偏移
		if (num_strand_correct_shapmers >= 2 * n)
		{
			if (strand == 1) // 相同的 strand
			{
				for (int i = 0; i <= num_strand_correct_shapmers - n - 1; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i + j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) <= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
						{
							score += -1;
						}
						else
						{
							score += 1;
						}
					}

					if (score > 0) // 认为这个 hapmer 是可信的，否则移除
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}

				// 处理最后几个 hapmers（这些 hapmers 没有 n 个后续 hapmers 来判断）
				for (int i = num_strand_correct_shapmers - n; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i - j];
						int qposi_p = hapmer_j.qposi;
						int tposi_p = hapmer_j.tposi;

						if ((qposi - qposi_p) <= 0 || (tposi - tposi_p) <= 0 || std::abs(std::abs(qposi - qposi_p) - std::abs(tposi - tposi_p)) > shift)
						{
							score += -1;
						}
						else
						{
							score += 1;
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
			//
			else // 不相同的 strand
			{
				for (int i = 0; i <= num_strand_correct_shapmers - n - 1; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i + j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) >= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
						{
							score += -1;
							// cout<<"Now score -1 : "<<score<<endl;
						}
						else
						{
							score += 1;
							// cout<<"Now score +1 : "<<score<<endl;
						}
					}

					if (score > 0) // 认为这个 hapmer 是可信的，否则移除
					{
						// cout<<"Accepted hapmer: "<<strand_correct_shapmers[i].hapmer<<", QPos: "<<strand_correct_shapmers[i].qposi<<", TPos: "<<strand_correct_shapmers[i].tposi<<", Strand: "<<strand_correct_shapmers[i].strand<<endl;
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}

				// 处理最后几个 hapmers（这些 hapmers 没有 n 个后续 hapmers 来判断）
				for (int i = num_strand_correct_shapmers - n; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 1; j <= n; ++j)
					{
						const auto &hapmer_j = strand_correct_shapmers[i - j];
						int qposi_p = hapmer_j.qposi;
						int tposi_p = hapmer_j.tposi;

						if ((qposi - qposi_p) <= 0 || (tposi - tposi_p) >= 0 || std::abs(std::abs(qposi - qposi_p) - std::abs(tposi - tposi_p)) > shift)
						{
							score += -1;
						}
						else
						{
							score += 1;
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
		}
		else if (num_strand_correct_shapmers >= 2 && num_strand_correct_shapmers < 2 * n)
		{
			if (strand == 1) // 相同的 strand
			{
				for (int i = 0; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 0; j < num_strand_correct_shapmers; ++j)
					{
						if (i == j)
						{
							continue;
						}

						const auto &hapmer_j = strand_correct_shapmers[j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if (j > i)
						{
							if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) <= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
						else
						{
							if ((qposi_n - qposi) >= 0 || (tposi_n - tposi) >= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
			else // 不同的 strand
			{
				for (int i = 0; i < num_strand_correct_shapmers; ++i)
				{
					const auto &hapmer_i = strand_correct_shapmers[i];
					int qposi = hapmer_i.qposi;
					int tposi = hapmer_i.tposi;

					int score = 0;
					for (int j = 0; j < num_strand_correct_shapmers; ++j)
					{
						if (i == j)
						{
							continue;
						}

						const auto &hapmer_j = strand_correct_shapmers[j];
						int qposi_n = hapmer_j.qposi;
						int tposi_n = hapmer_j.tposi;

						if (j > i)
						{
							if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) >= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
						else
						{
							if ((qposi_n - qposi) <= 0 || (tposi_n - tposi) <= 0 || std::abs(std::abs(qposi_n - qposi) - std::abs(tposi_n - tposi)) > shift)
							{
								score += -1;
							}
							else
							{
								score += 1;
							}
						}
					}

					if (score > 0)
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}
			}
		}
		else
		{
			continue;
		}
		if (posi_correct_shapmers.size() < min_num_shapmers)
		{
			continue;
		}

		int num_nonovlp_shapmers = 1; // 初始值为 1，考虑第一个 hapmer

		std::unordered_map<uint64_t, int> count;
		for (const auto &h : posi_correct_shapmers)
		{
			count[h.hapmer]++;
		}

		std::vector<Tread2shapmers> filtered;
		for (const auto &h : posi_correct_shapmers)
		{
			if (count[h.hapmer] == 1)
			{ // 只保留出现 1 次的
				filtered.push_back(h);
			}
		}
		posi_correct_shapmers.swap(filtered);
		if (filtered.size() < min_num_shapmers)
		{
			continue;
		}
		cout << "Quary name:" << s->name << endl;
		for (size_t i = 0; i < filtered.size(); ++i)
		{
			std::cout << "After filtering, Hapmer: " << filtered[i].hapmer
					  << "	, QPos: " << filtered[i].qposi
					  << "	, TPos: " << filtered[i].tposi
					  << "	, Strand: " << filtered[i].strand << "\n";
		}

		//...............................................................................

		for (size_t i = 0; i < filtered.size(); ++i)
		{
			const Tread2shapmers &hapmer1 = filtered[i];
			const Tread2shapmers &hapmer2 = filtered[i + 1];

			int qposi1 = hapmer1.qposi;
			int tposi1 = hapmer1.tposi;
			int qposi2 = hapmer2.qposi;
			int tposi2 = hapmer2.tposi;
			// cout << "  Filtered Hapmer: " << hapmer1.hapmer
			// 	 << ", QPos: " << hapmer1.qposi
			// 	 << ", TPos: " << hapmer1.tposi
			// 	 << ", Strand: " << hapmer1.strand << "\n";

			// 检查相邻 hapmers 是否满足距离要求
			if (std::abs(qposi2 - qposi1) >= aux->ch->k && std::abs(tposi2 - tposi1) >= aux->ch->k)
			{
				num_nonovlp_shapmers++;
			}
		}

		if (posi_correct_shapmers.size() > 60)
		{
			int num = posi_correct_shapmers.size() / 60;
			num_nonovlp_shapmers += num;
		}

		const Tread2shapmers &first_hapmer = filtered.front();
		const Tread2shapmers &last_hapmer = filtered.back();
		int qposi1 = first_hapmer.qposi;
		int tposi1 = first_hapmer.tposi;
		int qposi2 = last_hapmer.qposi;
		int tposi2 = last_hapmer.tposi;
		strand = first_hapmer.strand;
		int qstart, qend, tstart, tend;

		if (strand == 1)
		{ // same strand
			qstart = qposi1;
			qend = qposi2;
			tstart = tposi1;
			tend = tposi2;
		}
		else
		{ // different strand
			qstart = qposi1;
			qend = qposi2;
			tstart = tposi2;
			tend = tposi1;
		}
		// if (min(qend-qstart, tend-tstart)/max(qend-qstart, tend-tstart)) < 0.5:
		//     continue
		// // // 检查重叠比例
		double overlap = static_cast<double>(std::min(qend - qstart, tend - tstart)) / std::max(qend - qstart, tend - tstart);
		// cout << "	overlap ratio: " << overlap << " " << static_cast<double>(std::min(qend - qstart, tend - tstart)) << " " << std::max(qend - qstart, tend - tstart) << endl;
		if (num_nonovlp_shapmers < 2)
		{
			continue;
		}

		if (static_cast<double>(std::min(qend - qstart, tend - tstart)) / std::max(qend - qstart, tend - tstart) < 0.5)
		{
			continue;
		}

		// cout << "Num non-overlapping shapmers: " << num_nonovlp_shapmers << "\n";
		counts[target_id] = num_nonovlp_shapmers;
		forward[target_id] = strand;

		std::ofstream outFile;
		outFile.open("output.csv", std::ios::app);
		outFile << s->name << ","
				<< s->l_seq << ","
				<< qstart << ","
				<< qend << ","
				<< strand << ","
				<< target_id << ","
				<< tstart << ","
				<< tend << ","
				<< num_nonovlp_shapmers << ","
				<< filtered.size() << ","
				<< posi_correct_shapmers.size() << "\n";

		// 关闭文件
		outFile.close();

		int tlen = 0;
		// auto it = g_read_length_map.find(target_id);
		// if (it != g_read_length_map.end())
		// 	tlen = it->second;

		int mapping_quality = 0;
		if (!filtered.empty())
			mapping_quality = std::min(60, int(60.0 * num_nonovlp_shapmers / filtered.size()));

		// // ✅ 输出 PAF 格式
		// std::ofstream outFile_paf("output.paf", std::ios::app);
		// outFile_paf << s->name << "\t"				// query_name
		// 		<< s->l_seq << "\t"				// query_length
		// 		<< qstart << "\t"				// query_start
		// 		<< qend << "\t"					// query_end
		// 		<< strand << "\t"				// strand
		// 		<< target_id << "\t"			// target_name
		// 		<< tlen << "\t"					// target_length
		// 		<< tstart << "\t"				// target_start
		// 		<< tend << "\t"					// target_end
		// 		<< num_nonovlp_shapmers << "\t" // num_matching_bases
		// 		<< filtered.size() << "\t"		// alignment_block_length
		// 		<< mapping_quality<< "\n";			// mapping_quality
		// outFile_paf.close();
		// // （可选）附加标签
		// outFile << "\tov:i:" << overlap						 // 自定义字段 overlap
		// 		<< "\tpc:i:" << posi_correct_shapmers.size() // 正确 shapmers 数
		// 		<< std::endl;
	}

	uint16_t max_idx = -1;
	int max_val = 0;
	for (auto idx : counts)
	{
		if (max_val < idx.second)
		{
			max_val = idx.second;
			max_idx = idx.first;
		}
	}
}

static void *hifi_pipeline(void *shared, int step, void *_data)
{

	tb_shared_hifi_t *aux = (tb_shared_hifi_t *)shared;

	if (step == 0)
	{
		tb_step_hifi *s;
		s = (tb_step_hifi *)calloc(1, sizeof(tb_step_hifi));
		s->seq = bseq_read(aux->fp, CHUNK_SIZE, 0, &s->n_seq);
		s->aux = aux;
		if (s->n_seq)
		{
			s->mappings = (uint16_t *)calloc(s->n_seq, sizeof(uint16_t));
			s->mappings_forward = (bool *)calloc(s->n_seq, sizeof(bool));

			fprintf(stderr, "[M::%s] reads %d sequences\n", __func__, s->n_seq);
			return s;
		}
		else
			free(s);
	}
	else if (step == 1)
	{

		int i;
		tb_step_hifi *s = (tb_step_hifi *)_data;
		int success_mapping_num = 0;
		kt_for(aux->n_threads, hifi_worker, s, s->n_seq);
		// for (i = 0; i < s->n_seq; ++i)
		// {
		// 	mapping_res_t res;
		// 	// res.name = strdup(s->seq[i].name);
		// 	res.unit_id = s->mappings[i];
		// 	res.forward = s->mappings_forward[i];
		// 	if (s->mappings[i] != 65535)
		// 	{
		// 		success_mapping_num++;
		// 	}
		// 	aux->mappings->push_back(res);
		// 	free(s->seq[i].name);
		// 	free(s->seq[i].seq);
		// 	free(s->seq[i].qual);
		// 	free(s->seq[i].comment);
		// }
		// cout << "success_mapping_num: " << success_mapping_num << endl;

		double success_ratio = success_mapping_num / s->n_seq;
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences;\n", __func__,
				yak_realtime(), yak_cputime() / yak_realtime(), s->n_seq);

		if (s->seq)
			free(s->seq);
		if (s->mappings)
			free(s->mappings);
		if (s)
			free(s);
	}
	return 0;
}

static void *porec_pipeline(void *shared, int step, void *_data)
{

	tb_shared_t *aux = (tb_shared_t *)shared;

	if (step == 0)
	{
		tb_step_t *s;
		s = (tb_step_t *)calloc(1, sizeof(tb_step_t));
		s->seq = bseq_read(aux->fp, CHUNK_SIZE, 0, &s->n_seq);
		s->aux = aux;
		if (s->n_seq)
		{
			s->mappings = (uint16_t *)calloc(s->n_seq, sizeof(uint16_t));
			s->mappings_forward = (bool *)calloc(s->n_seq, sizeof(bool));

			fprintf(stderr, "[M::%s] contig %d sequences\n", __func__, s->n_seq);
			return s;
		}
		else
			free(s);
	}
	else if (step == 1)
	{

		int i;
		tb_step_t *s = (tb_step_t *)_data;
		int success_mapping_num = 0;
		kt_for(aux->n_threads, porec_worker, s, s->n_seq);
		for (i = 0; i < s->n_seq; ++i)
		{
			mapping_res_t res;
			// res.name = strdup(s->seq[i].name);
			res.unit_id = s->mappings[i];
			res.forward = s->mappings_forward[i];
			if (s->mappings[i] != 65535)
			{
				success_mapping_num++;
			}
			aux->mappings->push_back(res);
			free(s->seq[i].name);
			free(s->seq[i].seq);
			free(s->seq[i].qual);
			free(s->seq[i].comment);
		}
		cout << "success_mapping_num: " << success_mapping_num << endl;

		double success_ratio = success_mapping_num / s->n_seq;
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences;success_ratio %.2f .\n", __func__,
				yak_realtime(), yak_cputime() / yak_realtime(), s->n_seq, success_ratio);

		if (s->seq)
			free(s->seq);
		if (s->mappings)
			free(s->mappings);
		if (s)
			free(s);
	}
	return 0;
}

void hifi_do_mapping(pldat_t *pl1, pldat_t *pl2, bseq_file_t *hifi_fn, char *out_fn)
{
	if (hifi_fn == 0)
	{
		fprintf(stderr, "ERROR: Please give hifi files\n");
		exit(1);
	}

	ketopt_t o = KETOPT_INIT;
	int i, c, min_cnt = 2, mid_cnt = 5;
	tb_shared_hifi_t aux;

	memset(&aux, 0, sizeof(tb_shared_hifi_t));

	aux.n_threads = pl1->opt->n_thread, aux.print_diff = 0;
	aux.ratio_thres = 0.33;
	aux.k = pl1->h->k;
	aux.ch1 = pl1->h;
	aux.ch2 = pl2->h;

	aux.h_pos1 = pl1->h_pos; // 新增
	aux.h_pos2 = pl2->h_pos; // 新增
	aux.contigs_names_hap1 = pl1->names;
	aux.contigs_names_hap2 = pl2->names;
	aux.record_num = pl1->global_counter > pl2->global_counter ? pl1->global_counter : pl2->global_counter;

	cout << "Start mapping, total largerast hap number: " << aux.record_num << endl;
	aux.fp = hifi_fn;
	std::vector<mapping_res_t> *map1 = new std::vector<mapping_res_t>();

	aux.mappings = map1;

	aux.buf = (tb_buf_t *)calloc(aux.n_threads, sizeof(tb_buf_t));
	cout << "Start mapping " << endl;

	//....................................
	std::ofstream outFile;
	outFile.open("output.csv", std::ios::app);
	outFile << "name,l_seq,qstart,qend,forward,target_id,"
			   "tstart,tend,num_nonovlp_shapmers,"
			   "posi_correct_shapmers,num_shapmers\n";
	outFile.close();
	//......................................
	kt_pipeline(2, hifi_pipeline, &aux, 2);
	bseq_close(aux.fp);
	cout << "End hic1 mapping" << endl;
}

void porec_do_mapping(pldat_t *pl, bseq_file_t *hic_fn1, bseq_file_t *hic_fn2, char *out_fn)
{
	if (hic_fn1 == 0 || hic_fn2 == 0)
	{
		fprintf(stderr, "ERROR: Please give two hic files\n");
		exit(1);
	}

	ketopt_t o = KETOPT_INIT;
	int i, c, min_cnt = 2, mid_cnt = 5;
	tb_shared_t aux;

	memset(&aux, 0, sizeof(tb_shared_t));
	aux.n_threads = pl->opt->n_thread, aux.print_diff = 0;
	aux.ratio_thres = 0.33;
	aux.k = pl->h->k;
	aux.ch = pl->h;
	aux.h_pos = pl->h_pos; // 新增
	aux.record_num = pl->global_counter;
	cout << "Start mapping, total unitig number: " << aux.record_num << endl;
	aux.fp = hic_fn1;
	std::vector<mapping_res_t> *map1 = new std::vector<mapping_res_t>();
	aux.mappings = map1;
	aux.buf = (tb_buf_t *)calloc(aux.n_threads, sizeof(tb_buf_t));
	cout << "Start hic1 mapping " << endl;

	//....................................
	std::ofstream outFile;
	outFile.open("output.csv", std::ios::app);
	outFile << "name,l_seq,qstart,qend,forward,target_id,"
			   "tstart,tend,num_nonovlp_shapmers,"
			   "posi_correct_shapmers,num_shapmers\n";
	outFile.close();
	//......................................
	kt_pipeline(2, porec_pipeline, &aux, 2);
	bseq_close(aux.fp);
	cout << "End hic1 mapping" << endl;
	exit(1);
	for (i = 0; i < aux.n_threads; ++i)
	{
		free(aux.buf[i].s);
	}
	free(aux.buf);

	aux.fp = hic_fn2;
	std::vector<mapping_res_t> *map2 = new std::vector<mapping_res_t>();
	aux.mappings = map2;
	aux.buf = (tb_buf_t *)calloc(aux.n_threads, sizeof(tb_buf_t));
	cout << "Start hic2 mapping " << endl;
	kt_pipeline(2, porec_pipeline, &aux, 2);
	bseq_close(aux.fp);
	for (i = 0; i < aux.n_threads; ++i)
	{
		free(aux.buf[i].s);
	}
	free(aux.buf);

	printf("Map1 size: %ld;\n", map1->size());
	printf("Map2 size: %ld;\n", map2->size());
}

static void worker_for_polishing_target(void *data, long i, int tid) // callback for kt_for()
{
	stepdat_t *s = (stepdat_t *)data;
	yak_ch_t *h = s->p->h;
	yak_ch_t *h_pos = s->p->h_pos; // 新增
	ch_buf_t *b = &s->buf[i];
	b->n_ins += yak_ch_insert_list_kmer_full_for_polishing_target(h, h_pos, s->p->create_new, b->n, b->a, b->r, b->pos, b->f, i);
}

static void worker_for_polishing(void *data, long i, int tid) // callback for kt_for()
{
	stepdat_t *s = (stepdat_t *)data;
	yak_ch_t *h = s->p->h;
	yak_ch_t *h_pos = s->p->h_pos; // 新增
	ch_buf_t *b = &s->buf[i];
	b->n_ins += yak_ch_insert_list_kmer_full_for_polishing(h, h_pos, s->p->create_new, b->n, b->a, b->r, b->pos, b->f, i);
}

static void count_seq_buf_polishing_target(ch_buf_t *buf, int k, int p, int len, const char *seq, uint64_t seq_id) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	bool f;
	uint64_t x[2], mask = (1ULL << k * 2) - 1, shift = (k - 1) * 2;
	// char *qbuf_f = (char *)malloc(sizeof(char) * k); // 正向 k-mer 的质量值

	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i)
	{
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4)
		{												   // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;				   // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift; // reverse strand
			// 保存质量值到正向和反向缓冲区
			// qbuf_f[l % k] = qual[i];
			std::string kmer;
			if (++l >= k)
			{ // we find a k-mer

				f = x[0] < x[1];
				uint64_t y = f ? x[0] : x[1];
				ch_insert_buf(buf, p, yak_hash64(y, mask), seq_id, f, l - k); // 新增参数
																			  // }
			}
		}
		else
		{
			l = 0;
			x[0] = x[1] = 0; // if there is an "N", restart
		}
	}

	// free(qbuf_f); // 释放正向质量值缓冲区
}

static void count_seq_buf_polishing(ch_buf_t *buf, int k, int p, int len, const char *seq, const char *qual, uint64_t seq_id) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	bool f;
	uint64_t x[2], mask = (1ULL << k * 2) - 1, shift = (k - 1) * 2;
	// char *qbuf_f = (char *)malloc(sizeof(char) * k); // 正向 k-mer 的质量值

	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i)
	{
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4)
		{												   // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;				   // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift; // reverse strand
			// 保存质量值到正向和反向缓冲区
			// qbuf_f[l % k] = qual[i];
			std::string kmer;
			if (++l >= k)
			{ // we find a k-mer

				f = x[0] < x[1];
				uint64_t y = f ? x[0] : x[1];
				if (qual[i] == 'H')
				{
					// cout << "Found kmer with H quality: " << decode_kmer(x[0], k) << " ,pos:" << l - k << "\n";
					//  if(seq_id > 65535)
					//  {
					//  	fprintf(stderr, "ERROR: seq_id %d exceed the limit 65535\n", seq_id);
					//  	// exit(1);
					//  }
					//  uint64_t y2= yak_hash64(y, mask);
					//  cout<<"y2: "<<y2<<endl;
					ch_insert_buf(buf, p, yak_hash64(y, mask), seq_id, f, l - k); // 新增参数
				}
			}
		}
		else
		{
			l = 0;
			x[0] = x[1] = 0; // if there is an "N", restart
		}
	}

	// free(qbuf_f); // 释放正向质量值缓冲区
}

static void *worker_pipeline_target(void *data, int step, void *in) // callback for kt_pipeline()
{
	pldat_t *p = (pldat_t *)data;
	if (step == 0)
	{ // step 1: read a block of sequences
		int ret;
		stepdat_t *s;
		CALLOC(s, 1);
		s->p = p;
		s->global_bias = p->global_counter;
		while ((ret = kseq_read(p->ks)) >= 0)
		{
			p->global_counter++;
			// cout<<"Reading sequence number: "<<p->global_counter<<endl;
			int l = p->ks->seq.l;
			if (l < p->opt->k)
				continue;
			if (s->n == s->m)
			{
				s->m = s->m < 16 ? 16 : s->m + (s->n >> 1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			MALLOC(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			s->p->g_read_length_map[p->global_counter - 1] = l; // 记录读取的序列长度
			p->names->push_back(strdup(p->ks->name.s));
			cout << "Reading sequence name: " << p->ks->name.s << ", length: " << l << ",index :" << p->global_counter - 1 << endl;
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->opt->k + 1;

			// cout<<"Current Read Length: "<<l<<", Total K-mers: "<<s->nk<<endl;
			if (s->sum_len >= p->opt->chunk_size)
				break;
		}
		// cout<<"Total Read Sequences: "<<p->global_counter<<endl;
		if (s->sum_len == 0)
			free(s);
		else
			return s;
	}
	else if (step == 1)
	{ // step 2: extract k-mers
		stepdat_t *s = (stepdat_t *)in;
		int i, n = 1 << p->opt->pre, m;
		CALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i)
		{
			s->buf[i].m = m;
			MALLOC(s->buf[i].a, m);
			MALLOC(s->buf[i].r, m);
			MALLOC(s->buf[i].f, m);
			MALLOC(s->buf[i].pos, m); // 新增：分配pos数组
		}
		// cout<<"Start counting k-mers for polishing target..."<<endl;
		for (i = 0; i < s->n; ++i)
		{
			// cout<<"Counting k-mers for polishing target "<<i<<"..."<<endl;

			if (p->opt->k < 32)
				count_seq_buf_polishing_target(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i], s->global_bias + i);

			// cout<<"Finished counting k-mers for polishing target "<<i<<"..."<<endl;
			//  else
			//  	count_seq_buf_long(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i], s->global_bias + i);
			free(s->seq[i]);
		}
		free(s->seq);
		free(s->len);
		return s;
	}
	else if (step == 2)
	{ // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t *)in;
		int i, n = 1 << p->opt->pre;
		uint64_t n_ins = 0;
		// cout<<"Start inserting k-mers for polishing target into hash table..."<<endl;
		kt_for(p->opt->n_thread, worker_for_polishing_target, s, n);
		for (i = 0; i < n; i++)
		{
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
			free(s->buf[i].r);
			free(s->buf[i].f);
			free(s->buf[i].pos);
		}
		p->h->tot += n_ins;
		free(s->buf);
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences; %ld distinct k-mers in the hash table\n", __func__,
				yak_realtime(), yak_cputime() / yak_realtime(), s->n, (long)p->h->tot);
		fprintf(stderr, "Total Processed %ld sequences\n", p->global_counter);
		free(s);
	}
	return 0;
}
static void *worker_pipeline_hamper(void *data, int step, void *in) // callback for kt_pipeline()
{
	pldat_t *p = (pldat_t *)data;
	if (step == 0)
	{ // step 1: read a block of sequences
		int ret;
		stepdat_t *s;
		CALLOC(s, 1);
		s->p = p;
		s->global_bias = p->global_counter;
		s->qual = NULL; // 初始化质量值数组为 NULL
		while ((ret = kseq_read(p->ks)) >= 0)
		{
			p->global_counter++;
			int l = p->ks->seq.l;
			int ql = p->ks->qual.l; // 质量值的长度
			if (l < p->opt->k || l != ql)
				continue; // 确保长度一致
			if (l < p->opt->k)
				continue;
			if (s->n == s->m)
			{
				s->m = s->m < 16 ? 16 : s->m + (s->n >> 1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
				REALLOC(s->qual, s->m); // 重新分配质量值数组
			}
			MALLOC(s->seq[s->n], l);
			MALLOC(s->qual[s->n], ql); // 为质量值分配空间
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			memcpy(s->qual[s->n], p->ks->qual.s, ql); // 复制质量值

			p->names->push_back(strdup(p->ks->name.s));

			// uint64_t read_id = s->global_bias + s->n; // 当前read的全局ID
			// {
			// 	g_read_length_map[read_id] = l;
			// }

			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->opt->k + 1;
			if (s->sum_len >= p->opt->chunk_size)
				break;
		}
		if (s->sum_len == 0)
			free(s);
		else
			return s;
	}
	else if (step == 1)
	{ // step 2: extract k-mers
		stepdat_t *s = (stepdat_t *)in;
		int i, n = 1 << p->opt->pre, m;
		CALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i)
		{
			s->buf[i].m = m;
			MALLOC(s->buf[i].a, m);
			MALLOC(s->buf[i].r, m);
			MALLOC(s->buf[i].f, m);
			MALLOC(s->buf[i].pos, m); // 新增：分配pos数组
		}
		for (i = 0; i < s->n; ++i)
		{
			if (p->opt->k < 32)
				count_seq_buf_polishing(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i], s->qual[i], s->global_bias + i);
			// else
			// 	count_seq_buf_long(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i], s->global_bias + i);
			free(s->seq[i]);
			free(s->qual[i]); // 释放质量值内存
		}
		free(s->seq);
		free(s->qual); // 释放质量值数组内存
		free(s->len);
		return s;
	}
	else if (step == 2)
	{ // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t *)in;
		int i, n = 1 << p->opt->pre;
		uint64_t n_ins = 0;

		kt_for(p->opt->n_thread, worker_for_polishing, s, n);
		for (i = 0; i < n; i++)
		{
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
			free(s->buf[i].r);
			free(s->buf[i].f);
			free(s->buf[i].pos);
		}
		p->h->tot += n_ins;
		free(s->buf);
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %d sequences; %ld distinct k-mers in the hash table\n", __func__,
				yak_realtime(), yak_cputime() / yak_realtime(), s->n, (long)p->h->tot);
		fprintf(stderr, "Total Processed %ld sequences\n", p->global_counter);
		free(s);
	}
	return 0;
}
pldat_t *hamper_new1(const char *fn, const yak_copt_t *opt, yak_ch_t *h0)
{
	pldat_t *pl;
	CALLOC(pl, 1);
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0)
		return 0;
	pl->ks = kseq_init(fp);
	pl->opt = opt;
	if (h0)
	{
		pl->h = h0, pl->create_new = 0;
		assert(h0->k == opt->k && h0->pre == opt->pre);
	}
	else
	{
		pl->create_new = 1;
		pl->h = yak_ch_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	}
	pl->names = new std::vector<char *>();
	pl->h_pos = yak_ch_init(opt->k, 30, opt->bf_n_hash, opt->bf_shift); // 新增：初始化存pos的哈希表
	cout << "Initialize h_pos" << endl;
	kt_pipeline(3, worker_pipeline_hamper, pl, 3);
	kseq_destroy(pl->ks);
	gzclose(fp);
	return pl;
}

pldat_t *target_h1(const char *fn, const yak_copt_t *opt, yak_ch_t *h0)
{
	pldat_t *pl = new pldat_t();

	// pldat_t *pl;
	// CALLOC(pl, 1);
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0)
		return 0;
	pl->ks = kseq_init(fp);
	pl->opt = opt;
	if (h0)
	{
		pl->h = h0, pl->create_new = 0;
		assert(h0->k == opt->k && h0->pre == opt->pre);
	}
	else
	{
		pl->create_new = 1;
		pl->h = yak_ch_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	}
	pl->names = new std::vector<char *>();
	pl->h_pos = yak_ch_init(opt->k, 30, opt->bf_n_hash, opt->bf_shift); // 新增：初始化存pos的哈希表
	cout << "Initialize h_pos" << endl;
	if (!pl->ks)
	{
		fprintf(stderr, "ERROR: failed to init kseq from %s\n", fn);
		return 0;
	}

	kt_pipeline(3, worker_pipeline_target, pl, 3);
	kseq_destroy(pl->ks);
	gzclose(fp);
	return pl;
}

// pldat_t *pl;
// 	CALLOC(pl, 1);
// 	gzFile fp;
// 	if ((fp = gzopen(fn, "r")) == 0)
// 		return 0;
// 	pl->ks = kseq_init(fp);
// 	pl->opt = opt;
// 	if (h0)
// 	{
// 		pl->h = h0, pl->create_new = 0;
// 		assert(h0->k == opt->k && h0->pre == opt->pre);
// 	}
// 	else
// 	{
// 		pl->create_new = 1;
// 		pl->h = yak_ch_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
// 	}
// 	pl->h_pos = yak_ch_init(opt->k, 30, opt->bf_n_hash, opt->bf_shift); // 新增：初始化存pos的哈希表
// 	pl->names = new std::vector<char *>();

// 	kt_pipeline(3, worker_pipeline, pl, 3);
// 	kseq_destroy(pl->ks);
// 	gzclose(fp);

// 	return pl;

int main_polishing(int argc, char *argv[])
{
	pldat_t *h1, *h2; // 存储Hi-C数据的处理状态及信息
	int c;			  // 处理命令行参数
	char *fn_out = 0; // 输出文件名

	yak_copt_t opt; // 存储命令行选项和参数的信息。这些结构体用于配置Hi-C数据的映射过程。
	ketopt_t o = KETOPT_INIT;
	yak_copt_init(&opt); // 初始化结构体opt
	opt.pre = YAK_COUNTER_BITS;
	cout << "YAK_COUNTER_BITS: " << YAK_COUNTER_BITS << endl;
	opt.n_thread = 32;
	while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:o:m:", 0)) >= 0)
	{ // 解析后续所有参数
		if (c == 'k')
			opt.k = atoi(o.arg);
		else if (c == 'p')
			opt.pre = atoi(o.arg);
		else if (c == 'K')
			opt.chunk_size = atoi(o.arg);
		else if (c == 't')
			opt.n_thread = atoi(o.arg);
		else if (c == 'b')
			opt.bf_shift = atoi(o.arg);
		else if (c == 'o')
			fn_out = o.arg;
		else if (c == 'm')
			opt.mapping = atoi(o.arg);
	}
	if (argc - o.ind < 1)
	{ // 设置结构体opt变量，展示
		// Test
		fprintf(stderr, "main_polishing_test argc: %d\n", argc);
		fprintf(stderr, "o.ind: %d\n", o.ind);
		fprintf(stderr, " %s\n", argv[o.ind]);
		fprintf(stderr, "Using porec.cpp\n");
		fprintf(stderr, "Usage: pstools main_polishing [options] target_contigs1 target_contigs2 classpro_reads\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", opt.pre);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", opt.bf_shift);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "  -m INT     use haphic\n");
		fprintf(stderr, "  -o FILE    save mapping relationship to FILE []\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		return 1;
	}
	if (opt.pre < YAK_COUNTER_BITS)
	{
		fprintf(stderr, "ERROR: -p should be at least %d\n", YAK_COUNTER_BITS);
		return 1;
	}
	if (opt.k >= 64)
	{
		fprintf(stderr, "ERROR: -k must be smaller than 64\n");
		return 1;
	}
	else if (opt.k >= 32)
	{
		fprintf(stderr, "WARNING: counts are inexact if -k is greater than 31\n");
	}
	char *classpro_reads = argv[o.ind + 2]; // classpro reads

	// // 读取 GFA 文件
	// asg_t *graph = gfa_read(gfa_filename);
	// if (graph == NULL) {
	//     fprintf(stderr, "ERROR: failed to read GFA file: %s\n", gfa_filename);
	//     return 1;
	// }
	// 打开 reads 文件
	// bseq_file_t *classpro_file = bseq_open(classpro_fn);
	char *target_reads1 = argv[o.ind];	   // query reads
	char *target_reads2 = argv[o.ind + 1]; // query reads
	bseq_file_t *classpro_reads_fn = bseq_open(classpro_reads);

	cout << "main_polishing" << endl;
	// cout << "opt.pre : " << opt.pre << endl;
	h1 = target_h1(argv[o.ind], &opt, 0);
	g_read_length_map_hap1 = h1->g_read_length_map;
	h1->g_read_length_map.clear();

	h2 = target_h1(argv[o.ind + 1], &opt, 0);
	g_read_length_map_hap2 = h2->g_read_length_map;
	h2->g_read_length_map.clear();
	std::cout << "Start target_read" << std::endl;
	hifi_do_mapping(h1, h2, classpro_reads_fn, fn_out);
	std::cout << "End target_read" << std::endl;

	free(h1);
	free(h2);
	return 0;
}

int main_polishing_test(int argc, char *argv[])
{
	pldat_t *h;		  // 存储Hi-C数据的处理状态及信息
	int c;			  // 处理命令行参数
	char *fn_out = 0; // 输出文件名

	yak_copt_t opt; // 存储命令行选项和参数的信息。这些结构体用于配置Hi-C数据的映射过程。
	ketopt_t o = KETOPT_INIT;
	yak_copt_init(&opt); // 初始化结构体opt
	opt.pre = YAK_COUNTER_BITS;
	cout << "YAK_COUNTER_BITS: " << YAK_COUNTER_BITS << endl;
	opt.n_thread = 32;
	while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:o:m:", 0)) >= 0)
	{ // 解析后续所有参数
		if (c == 'k')
			opt.k = atoi(o.arg);
		else if (c == 'p')
			opt.pre = atoi(o.arg);
		else if (c == 'K')
			opt.chunk_size = atoi(o.arg);
		else if (c == 't')
			opt.n_thread = atoi(o.arg);
		else if (c == 'b')
			opt.bf_shift = atoi(o.arg);
		else if (c == 'o')
			fn_out = o.arg;
		else if (c == 'm')
			opt.mapping = atoi(o.arg);
	}
	if (argc - o.ind < 1)
	{ // 设置结构体opt变量，展示
		// Test
		fprintf(stderr, "main_polishing_test argc: %d\n", argc);
		fprintf(stderr, "o.ind: %d\n", o.ind);
		fprintf(stderr, " %s\n", argv[o.ind]);
		fprintf(stderr, "Using porec.cpp\n");
		fprintf(stderr, "Usage: pstools main_poreC [options] classpro_file qreads_file treads_file\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", opt.pre);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", opt.bf_shift);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "  -m INT     use haphic\n");
		fprintf(stderr, "  -o FILE    save mapping relationship to FILE []\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		return 1;
	}
	if (opt.pre < YAK_COUNTER_BITS)
	{
		fprintf(stderr, "ERROR: -p should be at least %d\n", YAK_COUNTER_BITS);
		return 1;
	}
	if (opt.k >= 64)
	{
		fprintf(stderr, "ERROR: -k must be smaller than 64\n");
		return 1;
	}
	else if (opt.k >= 32)
	{
		fprintf(stderr, "WARNING: counts are inexact if -k is greater than 31\n");
	}
	char *ONT_ref1_contigs = argv[o.ind + 1]; // Hi-C reads 1
	char *ONT_ref2_contigs = argv[o.ind + 2]; // Hi-C reads 2

	// // 读取 GFA 文件
	// asg_t *graph = gfa_read(gfa_filename);
	// if (graph == NULL) {
	//     fprintf(stderr, "ERROR: failed to read GFA file: %s\n", gfa_filename);
	//     return 1;
	// }
	// 打开 reads 文件
	// bseq_file_t *classpro_file = bseq_open(classpro_fn);
	char *query_fn = argv[o.ind]; // query reads
	bseq_file_t *hic_fn1 = bseq_open(ONT_ref1_contigs);
	bseq_file_t *hic_fn2 = bseq_open(ONT_ref2_contigs);

	cout << "main_polishing_test" << endl;
	// cout << "opt.pre : " << opt.pre << endl;
	h = hamper_new1(argv[o.ind], &opt, 0);
	std::cout << "Start target_read" << std::endl;

	std::cout << "End target_read" << std::endl;

	porec_do_mapping(h, hic_fn1, hic_fn2, fn_out);

	free(h);
	return 0;
}

int main_poreC_map_test(int argc, char *argv[])
{
	pldat_t *h;		  // 存储Hi-C数据的处理状态及信息
	int c;			  // 处理命令行参数
	char *fn_out = 0; // 输出文件名
	// char *gfa_filename = NULL; // GFA 文件名

	yak_copt_t opt; // 存储命令行选项和参数的信息
	ketopt_t o = KETOPT_INIT;
	yak_copt_init(&opt);
	opt.pre = YAK_COUNTER_BITS;
	opt.n_thread = 32;
	char *csv_filename = NULL;

	// =======================
	// 解析命令行参数
	// =======================
	while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:o:m:L:c:", 0)) >= 0)
	{
		if (c == 'k')
			opt.k = atoi(o.arg);
		else if (c == 'p')
			opt.pre = atoi(o.arg);
		else if (c == 'K')
			opt.chunk_size = atoi(o.arg);
		else if (c == 't')
			opt.n_thread = atoi(o.arg);
		else if (c == 'b')
			opt.bf_shift = atoi(o.arg);
		else if (c == 'o')
			fn_out = o.arg;
		else if (c == 'm')
			opt.mapping = atoi(o.arg);
		else if (c == 'L')
			opt.l = atoi(o.arg);
		else if (c == 'c')
			csv_filename = o.arg;
		// else if (c == 'g')
		//     gfa_filename = o.arg;
	}

	// =======================
	// 参数检查
	// =======================
	if (argc - o.ind < 3)
	{
		fprintf(stderr, "Usage: HapFold hic_mapping [options] utg.fa hic_reads_1 hic_reads_2 -o gfa_hic.out \n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -p INT     prefix length [%d]\n", opt.pre);
		fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", opt.bf_shift);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -K INT     chunk size [100m]\n");
		fprintf(stderr, "  -m INT     use haphic mapping\n");
		fprintf(stderr, "  -c FILE    input node_type CSV file\n");
		fprintf(stderr, "  -L INT     use model for phasing\n");
		fprintf(stderr, "  -o FILE    save mapping relationship to FILE []\n");
		fprintf(stderr, "  -g FILE    input GFA file (for debug)\n");
		fprintf(stderr, "Note: -b37 is recommended for human reads\n");
		return 1;
	}

	if (opt.pre < YAK_COUNTER_BITS)
	{
		fprintf(stderr, "ERROR: -p should be at least %d\n", YAK_COUNTER_BITS);
		return 1;
	}
	if (opt.k >= 64)
	{
		fprintf(stderr, "ERROR: -k must be smaller than 64\n");
		return 1;
	}
	else if (opt.k >= 32)
	{
		fprintf(stderr, "WARNING: counts are inexact if -k is greater than 31\n");
	}
	cerr<<" in main_poreC_map_test"<<endl;

	// =======================
	// 输入文件解析
	// =======================
	// char *classpro_fn = argv[o.ind];	  // class 文件
	char *query_fn = argv[o.ind];		  // query reads
	char *hic_fn1_name = argv[o.ind + 1]; // Hi-C reads 1
	char *hic_fn2_name = argv[o.ind + 2]; // Hi-C reads 2

	// // 读取 GFA 文件
	// asg_t *graph = gfa_read(gfa_filename);
	// if (graph == NULL) {
	//     fprintf(stderr, "ERROR: failed to read GFA file: %s\n", gfa_filename);
	//     return 1;
	// }
	// 打开 reads 文件
	// bseq_file_t *classpro_file = bseq_open(classpro_fn);
	bseq_file_t *query_read = bseq_open(query_fn);
	bseq_file_t *hic_fn1 = bseq_open(hic_fn1_name);
	bseq_file_t *hic_fn2 = bseq_open(hic_fn2_name);

	if (query_read == NULL || hic_fn1 == NULL || hic_fn2 == NULL)
	{
		fprintf(stderr, "ERROR: failed to open input read files\n");
		return 1;
	}

	if (csv_filename != NULL)
	{
		cout << "[INFO] Loading node types from " << csv_filename << " ..." << endl;
		load_node_type_csv(csv_filename, node_type_map1);
		std::cerr << "[INFO] Loaded " << node_type_map1.size()
				  << " node types from " << csv_filename << std::endl;
		cout << "[INFO] Loaded " << node_type_map1.size()
			 << " node types from " << csv_filename << std::endl;
	}

	// =======================
	// Hi-C 处理主逻辑
	// =======================
	h = yak_count_multi_new(argv[o.ind], &opt, 0);

	//.....
	do_mapping2(h, hic_fn1, hic_fn2, fn_out, node_type_map1);

	free(h->h);
	free(h);
	return 0;
}

