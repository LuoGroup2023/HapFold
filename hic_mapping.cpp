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
	uint32_t *pos;
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

typedef struct
{ // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	int create_new;
	kseq_t *ks;
	yak_ch_t *h;
	// yak_ch_t *h_pos; 
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
	// yak_ch_t *h_pos = s->p->h_pos; 
	ch_buf_t *b = &s->buf[i];

	// b->n_ins += yak_ch_insert_list_kmer_full(h, h_pos, s->p->create_new, b->n, b->a, b->r, b->pos, b->f, i);
	b->n_ins += yak_ch_insert_list_kmer_full(h, s->p->create_new, b->n, b->a, b->r, b->pos, b->f, i);
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
	// yak_ch_t *h_pos; 
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
	yak_ch_t *h_pos1;
	yak_ch_t *h_pos2;
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

std::string decode_kmer(uint64_t x, int k);

std::string decode_kmer(uint64_t x, int k)
{
	static const char nt4[] = {'A', 'C', 'G', 'T'};
	std::string kmer(k, 'A'); 
	for (int i = 0; i < k; ++i)
	{
		kmer[k - i - 1] = nt4[x & 0x3]; 
		x >>= 2;						
	}
	return kmer;
}
struct ReadState
{
	int last_q_pos = -1; 
	int last_t_pos = -1; 
	bool has_last = false;
	bool last_forward = false;
	int consist_count = 0; 
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
		REALLOC(b->pos, b->m);
	}
	b->a[b->n] = y;
	b->f[b->n] = f;
	// std::cout << "kmer: " << decode_kmer(y, 31) << " pos: " << kmer_pos << " seq_id: " << seq_id << std::endl;
	b->pos[b->n] = kmer_pos; 
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
				// uint64_t hashed_y = yak_hash64(y, mask);
				// static int debug_raw_cnt = 0;
				// // if (debug_raw_cnt < 20) {
				// std::cout << "[DEBUG-RAW] Real Kmer: " << decode_kmer(y, k)
				// 		  << " | Seq ID: " << seq_id
				// 		  << " | Pos: " << (i - k + 1)
				// 		  << " | f: " << f
				// 		  << " | Target Hash(a[j]): " << hashed_y << std::endl;
				// debug_raw_cnt++;
				// // }
				// // ==========================================
				ch_insert_buf(buf, p, yak_hash64(y, mask), seq_id, f, i - k + 1); 
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
				ch_insert_buf(buf, p, yak_hash64(y, mask), seq_id, f, i - k + 1); 
		}
		else
			l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
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
			MALLOC(s->buf[i].pos, m); 
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
	// pl->h_pos = yak_ch_init(opt->k, 30, opt->bf_n_hash, opt->bf_shift); 
	pl->names = new std::vector<char *>();

	kt_pipeline(3, worker_pipeline, pl, 3);
	kseq_destroy(pl->ks);
	gzclose(fp);

	return pl;
}

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
	// yak_ch_t *h_pos = aux->h_pos; 
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

				if (aux->ch->k < 32)
				{
					is_forward = x[0] < x[1];//forward is 1 , backward is 0
					y = yak_hash64(is_forward ? x[0] : x[1], mask);
					y2 = yak_hash64(is_forward ? x[1] : x[0], mask);
					// cout << "kmer=" << decode_kmer(x[0], aux->ch->k) << endl;
					// printf("[DEBUG-hic] hic -> HashKey: %lu \n",
					// 	   (unsigned long)y);
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
					uint16_t read_id = yak_ch_get_pos_1(aux->ch, y, &pos);
					// if (pos != 65535) 
					// {
					// 	static int debug_get_cnt = 0;
					// 	if (debug_get_cnt < 20)
					// 	{
					// 		printf("[DEBUG-GET] Queried k-mer: %s | Retrieved Contig ID: %u | Retrieved Pos: %u\n",
					// 			   decode_kmer(x[0], aux->ch->k).c_str(), (unsigned)read_id, (unsigned)pos);

					// 		__sync_fetch_and_add(&debug_get_cnt, 1);
					// 	}
					// }
					// // ==========================================
					// // cout<<"is_forward: "<<is_forward<<", stored_forward_bit: "<<((res & YAK_FORWARD_MASK) != 0)<<endl;
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
						hap.hapmer = x[0];											
						hap.qposi = l - aux->ch->k;									  
						hap.tposi = pos;											  
						hap.strand = !(((res & YAK_FORWARD_MASK) != 0) ^ is_forward);
						read2hapmers[target_id].push_back(hap);
						// printf("Found hapmer: %lu, QPos: %d, TPos: %d, Strand: %d, TargetID: %d\n",
						// 	   hap.hapmer, hap.qposi, hap.tposi, hap.strand, target_id);
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
		int num_same_strand = 0;
		int num_diff_strand = 0;

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
		int strand = (num_same_strand > num_diff_strand) ? 1 : 0;

		std::vector<Tread2shapmers> strand_correct_shapmers;
		for (const Tread2shapmers &shapmer : shapmers)
		{
			if (shapmer.strand == strand)
			{
				strand_correct_shapmers.push_back(shapmer);
			}
		}
		int num_strand_correct_shapmers = strand_correct_shapmers.size();
		if (num_strand_correct_shapmers < min_num_shapmers)
		{
			continue;
		}
		std::vector<Tread2shapmers> posi_correct_shapmers;
		int n = 3;		 
		int shift = 200; 
		if (num_strand_correct_shapmers >= 2 * n)
		{
			if (strand == 1) 
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

					if (score > 0) 
					{
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}

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
			else
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

					if (score > 0)
					{
						// cout<<"Accepted hapmer: "<<strand_correct_shapmers[i].hapmer<<", QPos: "<<strand_correct_shapmers[i].qposi<<", TPos: "<<strand_correct_shapmers[i].tposi<<", Strand: "<<strand_correct_shapmers[i].strand<<endl;
						posi_correct_shapmers.push_back(strand_correct_shapmers[i]);
					}
				}

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
			if (strand == 1) 
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
			else 
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
		int num_nonovlp_shapmers = 1; 
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
	// aux.h_pos = pl->h_pos; 
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

int main_poreC_map_test(int argc, char *argv[])
{
	pldat_t *h;		 
	int c;			 
	char *fn_out = 0;
	// char *gfa_filename = NULL; // GFA 文件名

	yak_copt_t opt; 
	ketopt_t o = KETOPT_INIT;
	yak_copt_init(&opt);
	opt.pre = YAK_COUNTER_BITS;
	opt.n_thread = 32;
	char *csv_filename = NULL;

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

	if (argc - o.ind < 3)
	{
		fprintf(stderr, "Usage: pstools porec_mapping [options] query_reads hic_reads_1 hic_reads_2 -g graph.gfa\n");
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
		fprintf(stderr, "  -g FILE    input GFA file (required)\n");
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
	cerr << " in main_poreC_map_test" << endl;
	// char *classpro_fn = argv[o.ind];	
	char *query_fn = argv[o.ind];		 
	char *hic_fn1_name = argv[o.ind + 1]; // Hi-C reads 1
	char *hic_fn2_name = argv[o.ind + 2]; // Hi-C reads 2

	// asg_t *graph = gfa_read(gfa_filename);
	// if (graph == NULL) {
	//     fprintf(stderr, "ERROR: failed to read GFA file: %s\n", gfa_filename);
	//     return 1;
	// }
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

	h = yak_count_multi_new(argv[o.ind], &opt, 0);
	long peak_bytes = yak_peakrss();
	double peak_mb = peak_bytes / 1024.0 / 1024.0;
	double peak_gb = peak_bytes / 1024.0 / 1024.0 / 1024.0;
	printf("\n[Memory Info]\n");
	printf("Peak RSS: %ld Bytes\n", peak_bytes);
	printf("Peak RSS: %.2f MB\n", peak_mb);
	printf("Peak RSS: %.2f GB\n", peak_gb);
	//.....
	do_mapping2(h, hic_fn1, hic_fn2, fn_out, node_type_map1);

	free(h->h);
	free(h);
	return 0;
}
