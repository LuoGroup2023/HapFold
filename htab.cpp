#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "kthread.h"
#include "yak-priv.h"
#include <vector>
#include <pthread.h>
#include <set>
#include <map>
#include "khashl.h"
#include <inttypes.h>
#include <iostream>
#include <mutex> // hash table
#include <unordered_map>
#define yak_ch_eq(a, b) ((a) >> YAK_COUNTER_BITS == (b) >> YAK_COUNTER_BITS) // lower 8 bits for counts; higher bits for k-mer
#define yak_ch_hash(a) ((a) >> YAK_COUNTER_BITS)
KHASHL_SET_INIT(, yak_ht_t, yak_ht, uint64_t, yak_ch_hash, yak_ch_eq)

std::unordered_map<uint64_t, std::vector<RepeatInfo>> repeated_kmers_map;

std::mutex repeated_mutex;

//........................TODO:

typedef struct
{
	uint32_t read_id;
	uint32_t pos;
	uint8_t strand; // 0/1
} yak_occ_t;

typedef struct
{
	uint64_t off; // offset into occ_pool
	uint32_t cnt; // number of entries
} yak_poslist_t;

// 使用 khashl 生成 kmer -> yak_poslist_t 的 map
// 我们用前缀 yak_pos （会生成 yak_pos_put, yak_pos_get, yak_pos_init, yak_pos_destroy 等）
KHASHL_MAP_INIT(KH_LOCAL, yak_ht_pos_t, yak_pos, uint64_t, yak_poslist_t, kh_hash_uint64, kh_eq_generic)
// 上面宏等价于：
// HType = yak_ht_pos_t
// prefix = yak_pos
// key type = uint64_t
// val type = yak_poslist_t
// hash fn = kh_hash_uint64（在 khashl.h 中已声明）
// eq fn = kh_eq_generic

typedef struct
{
	yak_ht_pos_t *h; // khashl map pointer (实际类型由宏定义)
	yak_occ_t *occ_pool;
	uint64_t occ_cap;
	uint64_t occ_size;
} yak_ch_pos_t;

// 初始化
yak_ch_pos_t *yak_ch_pos_init(size_t init_cap)
{
	yak_ch_pos_t *p = (yak_ch_pos_t *)calloc(1, sizeof(yak_ch_pos_t));
	if (!p)
		return NULL;
	p->h = yak_pos_init(); // 由宏生成的初始化函数
	p->occ_cap = init_cap ? init_cap : (1ULL << 20);
	p->occ_size = 0;
	p->occ_pool = (yak_occ_t *)malloc(sizeof(yak_occ_t) * p->occ_cap);
	if (!p->occ_pool)
	{
		free(p);
		return NULL;
	}
	return p;
}

// 销毁
void yak_ch_pos_destroy(yak_ch_pos_t *p)
{
	if (!p)
		return;
	if (p->occ_pool)
		free(p->occ_pool);
	if (p->h)
		yak_pos_destroy(p->h);
	free(p);
}

// 内部扩容 occ_pool
static inline void yak_ch_pos_ensure_cap(yak_ch_pos_t *p, uint64_t need)
{
	if (p->occ_size + need <= p->occ_cap)
		return;
	uint64_t newcap = p->occ_cap ? p->occ_cap : 1;
	while (p->occ_size + need > newcap)
		newcap <<= 1;
	p->occ_pool = (yak_occ_t *)realloc(p->occ_pool, sizeof(yak_occ_t) * newcap);
	p->occ_cap = newcap;
}

// 插入单个 occurrence 到表中
// kmer: 使用与 yak 中一致的 kmer 编码（应与主表 key 一致）
// read_id,pos,strand: occurrence 信息
static inline void
yak_ch_pos_insert_one(yak_ch_pos_t *p, uint64_t kmer,
					  uint32_t read_id, uint32_t pos, uint8_t strand)
{
	int absent;
	// 调用宏生成的 put （prefix_put）
	khint_t idx = yak_pos_put(p->h, kmer, &absent);
	yak_poslist_t *pl = &kh_val(p->h, idx); // kh_val 宏可用
	if (absent)
	{
		// 第一次见到 kmer，pl->off 要指向当前 pool 末尾
		pl->off = p->occ_size;
		pl->cnt = 0;
	}
	// 确保 pool 容量
	yak_ch_pos_ensure_cap(p, 1);
	// 写入 occurrence
	yak_occ_t *occ = &p->occ_pool[p->occ_size++];
	occ->read_id = read_id;
	occ->pos = pos;
	occ->strand = strand;
	pl->cnt++;
}

// 查询函数：返回 poslist（指向 pool 的片段）
// 返回 pl->cnt，并用 *out_ptr 设置为指向第一个 yak_occ_t（注意：返回的内存由 yak_ch_pos_t 管理）
static inline uint32_t yak_ch_pos_get(yak_ch_pos_t *p, uint64_t kmer, yak_occ_t **out_ptr)
{
	khint_t idx = yak_pos_get(p->h, kmer);
	if (idx == kh_end(p->h))
	{
		*out_ptr = NULL;
		return 0;
	}
	yak_poslist_t *pl = &kh_val(p->h, idx);
	*out_ptr = p->occ_pool + pl->off;
	return pl->cnt;
}

//........................................

yak_ch_t *yak_ch_init(int k, int pre, int n_hash, int n_shift)
{
	yak_ch_t *h;
	int i;
	if (pre < YAK_COUNTER_BITS)
		return 0;
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1 << h->pre);
	for (i = 0; i < 1 << h->pre; ++i)
		h->h[i].h = yak_ht_init();
	if (n_hash > 0 && n_shift > h->pre)
	{
		h->n_hash = n_hash, h->n_shift = n_shift;
		for (i = 0; i < 1 << h->pre; ++i)
			h->h[i].b = yak_bf_init(h->n_shift - h->pre, h->n_hash);
	}
	return h;
}

void yak_ch_destroy_bf(yak_ch_t *h)
{
	int i;
	for (i = 0; i < 1 << h->pre; ++i)
	{
		if (h->h[i].b)
			yak_bf_destroy(h->h[i].b);
		h->h[i].b = 0;
	}
}

void yak_ch_destroy(yak_ch_t *h)
{
	int i;
	if (h == 0)
		return;
	yak_ch_destroy_bf(h);
	for (i = 0; i < 1 << h->pre; ++i)
		yak_ht_destroy(h->h[i].h);
	free(h->h);
	free(h);
}

int yak_ch_insert_list(yak_ch_t *h, int create_new, int n, const uint64_t *a)
{
	int j, mask = (1 << h->pre) - 1, n_ins = 0;
	yak_ch1_t *g;
	if (n == 0)
		return 0;
	g = &h->h[a[0] & mask];
	for (j = 0; j < n; ++j)
	{
		int ins = 1, absent;
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		if ((a[j] & mask) != (a[0] & mask))
			continue;
		if (create_new)
		{
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);
			if (ins)
			{
				k = yak_ht_put(g->h, x << YAK_COUNTER_BITS, &absent);
				if (absent)
					++n_ins;
				if ((kh_key(g->h, k) & YAK_MAX_COUNT) < YAK_MAX_COUNT)
					++kh_key(g->h, k);
			}
		}
		else
		{
			k = yak_ht_get(g->h, x << YAK_COUNTER_BITS);
			if (k != kh_end(g->h) && (kh_key(g->h, k) & YAK_MAX_COUNT) < YAK_MAX_COUNT)
				++kh_key(g->h, k);
		}
	}
	return n_ins;
}

int yak_ch_insert_list_kmer_record_mapping(yak_ch_t *h, int create_new, int n, const uint64_t *a, const bool *f, recordset_t *recordset, const uint16_t *r, long i, std::set<uint64_t> *deleted)
{
	int j, mask = (1 << h->pre) - 1, n_ins = 0;
	yak_ch1_t *g;
	if (n == 0)
		return 0;
	g = &h->h[a[0] & mask];
	for (j = 0; j < n; ++j)
	{
		int ins = 1, absent;
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		if ((a[j] & mask) != (a[0] & mask))
			continue;
		if (create_new)
		{
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);
			if (ins)
			{
				k = yak_ht_put(g->h, x << YAK_COUNTER_BITS, &absent);
				if (absent)
				{
					++n_ins;
					kh_key(g->h, k) |= r[j];
					// printf("0x%016" PRIx64 "\n", x);
					if (f[j])
					{
						kh_key(g->h, k) |= YAK_FORWARD_MASK;
					}
				}
				else if (((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
				{
					kh_key(g->h, k) |= YAK_REPEAT_MASK;
				}
			}
		}
		else
		{
			k = yak_ht_get(g->h, x << YAK_COUNTER_BITS);
			if (k != kh_end(g->h) && ((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
			{
				kh_key(g->h, k) |= YAK_REPEAT_MASK;
			}
		}
	}
	return n_ins;
}

// int yak_ch_insert_list_kmer_record_mapping(yak_ch_t *h, int create_new, int n, const uint64_t *a, const bool *f, recordset_t *recordset, const uint16_t *r, long i, std::set<uint64_t> *deleted)
// {
// 	int j, mask = (1 << h->pre) - 1, n_ins = 0;
// 	yak_ch1_t *g;
// 	if (n == 0)
// 		return 0;
// 	g = &h->h[a[0] & mask];
// 	for (j = 0; j < n; ++j)
// 	{
// 		int ins = 1, absent;
// 		uint64_t x = a[j] >> h->pre;
// 		khint_t k;
// 		if ((a[j] & mask) != (a[0] & mask))
// 			continue;
// 		if (create_new)
// 		{
// 			if (g->b)
// 				ins = (yak_bf_insert(g->b, x) == h->n_hash);
// 			// printf("aj=%llu, x=%llu\n", a[j], x);
// 			// printf("ins: %d\n", ins);
// 			if (ins)
// 			{
// 				k = yak_ht_put(g->h, x << YAK_COUNTER_BITS, &absent);
// 				if (absent)
// 				{
// 					++n_ins;
// 					kh_key(g->h, k) |= r[j];
// 					//printf("insert readID=%u\n", r[j]);
// 					if (f[j])
// 					{
// 						kh_key(g->h, k) |= YAK_FORWARD_MASK;
// 					}
// 				}
// 				else if (((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
// 				{
// 					//printf("in else\n");
// 					kh_key(g->h, k) |= YAK_REPEAT_MASK;
// 				}
// 			}
// 		}
// 		else
// 		{
// 			k = yak_ht_get(g->h, x << YAK_COUNTER_BITS);
// 			if (k != kh_end(g->h) && ((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
// 			{
// 				kh_key(g->h, k) |= YAK_REPEAT_MASK;
// 			}
// 		}
// 	}

// 	return n_ins;
// }
extern std::string decode_kmer(uint64_t x, int k); // 声明decode_kmer

// int yak_hash_put(yak_ch_t *h)
// {

// 	for (khint_t i = 0; i < kh_end(h); i++)
// 	{
// 		if (kh_exist(h, i))
// 		{ // 检查该 bucket 是否有效
// 			printf("bucket[%u] key = 0x%lx\n", i, kh_key(h, i));
// 		}
// 	}
// }
int yak_ch_insert_list_kmer_full(yak_ch_t *h, yak_ch_t *h_pos,
								 int create_new, int n,
								 const uint64_t *a, const uint64_t *r,
								 const uint32_t *pos, const bool *f,
								 long i)
{
	int j, mask = (1 << h->pre) - 1, mask_pos = (1 << 30) - 1, n_ins = 0;
	yak_ch1_t *g, *g_pos;
	if (n == 0)
		return 0;
	// g = &h->h[a[0] & mask];
	// g_pos = &h_pos->h[a[0] & mask_pos];

	for (j = 0; j < n; ++j)
	{
		int ins = 1, absent;
		yak_ch1_t *g = &h->h[a[j] & mask];
		yak_ch1_t *g_pos = &h_pos->h[a[j] & mask_pos];
		uint64_t x = a[j] >> h->pre;
		uint64_t x2 = a[j] >> 30;
		khint_t k, k_pos;
		// std::cout<<"j: "<<j<<std::endl;
		// if ((a[j] & mask) != (a[0] & mask))
		// 	continue;
		// if ((a[j] & mask_pos) != (a[0] & mask_pos))
		// 	continue;
		if (create_new)
		{
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);

			if (ins)
			{
				// 插入 kmer 本体 + r[j] + forward mask
				k = yak_ht_put(g->h, x << YAK_COUNTER_BITS, &absent);
				if (absent)
				{
					// std::cout<<"g_pos : a[0] & mask_pos :" << (a[0] & mask_pos) << std::endl;
					// std::cout<<"g : a[0] & mask :" << (a[0] & mask) << std::endl;
					++n_ins;
					kh_key(g->h, k) |= r[j];
					if (f[j])
					{
						kh_key(g->h, k) |= YAK_FORWARD_MASK;
					}

					//....
					uint64_t key_h = kh_key(g->h, k);

					// int get_k = yak_ch_get_k(h, a[j]);
					// if (get_k == -1)
					// {
					// 	uint16_t utg = get_k & YAK_KEY_MASK;
					// 	// std::cout << "utg=" << utg << std::endl;
					// }
					// kmer 本体 (高位去掉计数位)
					uint64_t kmer = key_h >> YAK_COUNTER_BITS;
					std::string kmer_str_in_h = decode_kmer(kmer, k);
					// printf("raw_key=%lu, kmer_code=%lu, kmer=%s\n",
					// 	   key_h, kmer, kmer_str_in_h.c_str());
					// read id（在低位的 KEY_MASK 范围里）
					uint64_t read_id = key_h & YAK_KEY_MASK;

					// 是否 forward
					bool is_forward = key_h & YAK_FORWARD_MASK;

					// 是否重复
					bool is_repeat = key_h & YAK_REPEAT_MASK;
					std::string kmer_str = decode_kmer(a[j], h->k);
					// std::cout << a[j] << std::endl;
					// printf("pos=%u, readID=%lu, forward=%d, repeat=%d\n",
					// pos[j], read_id, is_forward, is_repeat);
					//...
					//  插入位置信息
					int absent1;
					k_pos = yak_ht_put(g_pos->h, x2 << 30, &absent1);
					// std::cout << absent1 << std::endl;
					if (absent1)
					{

						kh_key(g_pos->h, k_pos) |= (pos[j] & YAK_POS_MASK);
					}
					else
					{
						printf("absent1 = %d", absent1);
						int ha = (kh_key(g_pos->h, k_pos) & YAK_KEY_MASK) ^ r[j];
					}

					// std::cout << ((a[j] >> 30) << 30) << std::endl;
					k = yak_ht_get(g_pos->h, x2 << 30);
					// std::cout<<"x2 << 30: "<<(x2 << 30)<<std::endl;

					// std::cout<<"g_pos->count in: "<<g_pos->h->count<<" g_pos->keys: "<<g_pos->h->keys<<std::endl;
					//  std::cout<<k<<std::endl;
					//  std::cout<<k_pos<<std::endl;
					uint64_t key_pos = kh_key(g_pos->h, k);
					uint32_t position = key_pos & YAK_POS_MASK;
					uint32_t position2 = 65535;
					// printf("position=%u\n", position);
					// std::cout<<"a[j] & mask_pos: "<<(a[j] & mask_pos)<<std::endl;
					// std::cout<<"a[j] : "<<(a[j])<<std::endl;
					yak_ch_get_pos(h, h_pos, a[j], &position2);
					// printf("position2=%u\n", position2);
				}

				else if (((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
				{
					// printf("in else\n");
					kh_key(g->h, k) |= YAK_REPEAT_MASK;
					bool is_repeat = 1;
					// std::cout << "is_repeat: " << is_repeat << std::endl;
				}
			}
		}
		// } else {
		//     k = yak_ht_get(g->h, x << YAK_COUNTER_BITS);
		//     if (k != kh_end(g->h) && ((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
		//         kh_key(g->h, k) |= YAK_REPEAT_MASK;
		// }
	}
	return n_ins;
}

int yak_ch_insert_list_kmer_full_for_polishing_target(yak_ch_t *h, yak_ch_t *h_pos,
													  int create_new, int n,
													  const uint64_t *a, const uint64_t *r,
													  const uint32_t *pos, const bool *f,
													  long i)
{
	int j, mask = (1 << h->pre) - 1, mask_pos = (1 << 30) - 1, n_ins = 0;
	yak_ch1_t *g, *g_pos;
	if (n == 0)
		return 0;
	// g = &h->h[a[0] & mask];
	// g_pos = &h_pos->h[a[0] & mask_pos];

	for (j = 0; j < n; ++j)
	{
		int ins = 1, absent;
		yak_ch1_t *g = &h->h[a[j] & mask];
		yak_ch1_t *g_pos = &h_pos->h[a[j] & mask_pos];
		uint64_t x = a[j] >> h->pre;
		uint64_t x2 = a[j] >> 30;
		khint_t k, k_pos;
		// std::cout<<"j: "<<j<<std::endl;
		// if ((a[j] & mask) != (a[0] & mask))
		// 	continue;
		// if ((a[j] & mask_pos) != (a[0] & mask_pos))
		// 	continue;
		if (create_new)
		{
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);

			if (ins)
			{
				// 插入 kmer 本体 + r[j] + forward mask
				k = yak_ht_put(g->h, x << YAK_COUNTER_BITS, &absent);
				if (absent)
				{
					++n_ins;
					kh_key(g->h, k) |= r[j];
					if (f[j])
					{
						kh_key(g->h, k) |= YAK_FORWARD_MASK;
					}

					int absent1;
					k_pos = yak_ht_put(g_pos->h, x2 << 30, &absent1);
					// std::cout << absent1 << std::endl;
					if (absent1)
					{

						kh_key(g_pos->h, k_pos) |= (pos[j] & YAK_POS_MASK);
					}
					else
					{
						int ha = (kh_key(g_pos->h, k_pos) & YAK_KEY_MASK) ^ r[j];
					}
				}

				else if (((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
				{
					// printf("in else\n");
					kh_key(g->h, k) |= YAK_REPEAT_MASK;
				}
			}
		}
	}
	return n_ins;
}

int yak_ch_insert_list_kmer_full_for_polishing(yak_ch_t *h, yak_ch_t *h_pos,
											   int create_new, int n,
											   const uint64_t *a, const uint64_t *r,
											   const uint32_t *pos, const bool *f,
											   long i)
{
	int j, mask = (1 << h->pre) - 1, mask_pos = (1 << 30) - 1, n_ins = 0;
	yak_ch1_t *g, *g_pos;
	if (n == 0)
		return 0;
	// g = &h->h[a[0] & mask];
	// g_pos = &h_pos->h[a[0] & mask_pos];

	for (j = 0; j < n; ++j)
	{
		int ins = 1, absent;
		yak_ch1_t *g = &h->h[a[j] & mask];
		yak_ch1_t *g_pos = &h_pos->h[a[j] & mask_pos];
		uint64_t x = a[j] >> h->pre;
		uint64_t x2 = a[j] >> 30;
		khint_t k, k_pos;
		// std::cout<<"j: "<<j<<std::endl;
		// if ((a[j] & mask) != (a[0] & mask))
		// 	continue;
		// if ((a[j] & mask_pos) != (a[0] & mask_pos))
		// 	continue;
		if (create_new)
		{
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);

			if (ins)
			{
				// 插入 kmer 本体 + r[j] + forward mask
				k = yak_ht_put(g->h, x << YAK_COUNTER_BITS, &absent);
				if (absent)
				{
					// std::cout<<"g_pos : a[0] & mask_pos :" << (a[0] & mask_pos) << std::endl;
					// std::cout<<"g : a[0] & mask :" << (a[0] & mask) << std::endl;
					++n_ins;
					kh_key(g->h, k) |= r[j];
					if (f[j])
					{
						kh_key(g->h, k) |= YAK_FORWARD_MASK;
					}

					//....
					uint64_t key_h = kh_key(g->h, k);

					int get_k = yak_ch_get_k(h, a[j]);
					if (get_k != -1)
					{
						uint32_t utg = get_k & YAK_KEY_MASK;
						// std::cout << "utg=" << utg <<" a[j]: "<< a[j]<< std::endl;
						//  std::cout << a[j] << std::endl;
					}
					// kmer 本体 (高位去掉计数位)
					uint64_t kmer = key_h >> YAK_COUNTER_BITS;
					std::string kmer_str_in_h = decode_kmer(kmer, k);
					// printf("raw_key=%lu, kmer_code=%lu, kmer=%s\n",
					// 	   key_h, kmer, kmer_str_in_h.c_str());
					// read id（在低位的 KEY_MASK 范围里）
					uint64_t read_id = key_h & YAK_KEY_MASK;

					// 是否 forward
					bool is_forward = key_h & YAK_FORWARD_MASK;

					// 是否重复
					bool is_repeat = key_h & YAK_REPEAT_MASK;
					std::string kmer_str = decode_kmer(a[j], h->k);

					// printf("pos=%u, readID=%lu, forward=%d, repeat=%d\n",
					// pos[j], read_id, is_forward, is_repeat);
					//...
					//  插入位置信息
					int absent1;
					k_pos = yak_ht_put(g_pos->h, x2 << 30, &absent1);
					// std::cout << absent1 << std::endl;
					if (absent1)
					{

						kh_key(g_pos->h, k_pos) |= (pos[j] & YAK_POS_MASK);
					}
					else
					{
						printf("absent1 = %d", absent1);
						int ha = (kh_key(g_pos->h, k_pos) & YAK_KEY_MASK) ^ r[j];
					}

					// std::cout << ((a[j] >> 30) << 30) << std::endl;
					k = yak_ht_get(g_pos->h, x2 << 30);
					// std::cout<<"x2 << 30: "<<(x2 << 30)<<std::endl;

					// std::cout<<"g_pos->count in: "<<g_pos->h->count<<" g_pos->keys: "<<g_pos->h->keys<<std::endl;
					//  std::cout<<k<<std::endl;
					//  std::cout<<k_pos<<std::endl;
					uint64_t key_pos = kh_key(g_pos->h, k);
					uint32_t position = key_pos & YAK_POS_MASK;
					uint32_t position2 = 65535;
					// printf("position=%u\n", position);
					//  std::cout<<"a[j] & mask_pos: "<<(a[j] & mask_pos)<<std::endl;
					//  std::cout<<"a[j] : "<<(a[j])<<std::endl;
					yak_ch_get_pos(h, h_pos, a[j], &position2);

					// printf("position2=%u\n", position2);
				}

				else if (((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
				{
					// printf("in else\n");
					kh_key(g->h, k) |= YAK_REPEAT_MASK;
					bool is_repeat = 1;
					// std::cout << "is_repeat: " << is_repeat << std::endl;

					// 在并行循环里
					if (is_repeat)
					{
						// 保存额外信息到外部 map
						RepeatInfo info;
						info.read_name = r[j]; // 或 r[j] 如果你用 id
						info.pos = pos[j];
						info.forward = f[j];

						std::lock_guard<std::mutex> lock(repeated_mutex);
						repeated_kmers_map[a[j]].push_back(info);

						// std::cout << "Detected repeat kmer: " << a[j] << ", readID: " << r[j] << ", pos: " << pos[j] << ", forward: " << f[j] << "\n";
					}
				}
			}
		}
		// } else {
		//     k = yak_ht_get(g->h, x << YAK_COUNTER_BITS);
		//     if (k != kh_end(g->h) && ((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
		//         kh_key(g->h, k) |= YAK_REPEAT_MASK;
		// }
	}
	return n_ins;
}

int yak_ch_insert_list_kmer_record_mapping2(yak_ch_t *h, yak_ch_t *h_pos, int create_new, int n, const uint64_t *a, const bool *f, recordset_t *recordset, const uint16_t *r, const uint32_t *pos, long i, std::set<uint64_t> *deleted)
{
	int j, mask = (1 << h->pre) - 1, n_ins = 0;
	yak_ch1_t *g;
	if (n == 0)
		return 0;
	g = &h->h[a[0] & mask];
	for (j = 0; j < n; ++j)
	{
		int ins = 1, absent;
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		if ((a[j] & mask) != (a[0] & mask))
			continue;
		if (create_new)
		{
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);
			if (ins)
			{
				k = yak_ht_put(g->h, x << YAK_COUNTER_BITS, &absent);
				if (absent)
				{
					++n_ins;
					kh_key(g->h, k) |= r[j];
					if (f[j])
					{
						kh_key(g->h, k) |= YAK_FORWARD_MASK;
					}
					// ---- 新增：插入pos到h_pos ----
					yak_ch1_t *g_pos = &h_pos->h[a[0] & ((1 << h_pos->pre) - 1)];
					int absent_pos;
					khint_t k_pos = yak_ht_put(g_pos->h, a[j], &absent_pos);
					kh_key(g_pos->h, k_pos) = (pos[j] & YAK_POS_MASK); // 只存低30位
																	   // ---- 新增结束 ----
				}
				else if (((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
				{
					kh_key(g->h, k) |= YAK_REPEAT_MASK;
				}
			}
		}
		else
		{
			k = yak_ht_get(g->h, x << YAK_COUNTER_BITS);
			if (k != kh_end(g->h) && ((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
			{
				kh_key(g->h, k) |= YAK_REPEAT_MASK;
			}
		}
	}

	return n_ins;
}

// int yak_ch_insert_list_kmer_pos(yak_ch_t *h, yak_ch_t *h_pos, int create_new, int n, const uint64_t *a, const uint16_t *r, const uint32_t*pos, long i)
// {
//     int j, mask = (1<<h->pre) - 1, n_ins = 0;
//     yak_ch1_t *g;
//     yak_ch1_t *g_pos;
//     if (n == 0) return 0;
//     g = &h->h[a[0]&mask];
//     g_pos = &h_pos->h[a[0]&((1<<h_pos->pre)-1)];

//     for (j = 0; j < n; ++j) {
//         int ins = 1, absent;
// 		//std::cout<<"a: "<<a[j]<<" r[j]: "<<r[j]<<" pos[j]: "<<pos[j]<<std::endl;
// 		printf("a[j]: %llx\n", a[j]);
// 		printf("r[j]: %d\n", r[j]);
// 		printf("pos[j]: %u\n", pos[j]);
//         uint64_t x = a[j] >> h->pre;
//         khint_t k;
//         if ((a[j]&mask) != (a[0]&mask)) continue;
//         if (create_new) {
//             if (g->b)
//                 ins = (yak_bf_insert(g->b, x) == h->n_hash);
//             if (ins) {
//                 k = yak_ht_put(g->h, (a[j] >> h->pre)<<YAK_COUNTER_BITS, &absent);
//                 if (absent){
//                     kh_key(g->h, k)|=r[j];
//                     // 用完整的k-mer作为key，直接赋值
//                     k = yak_ht_put(g_pos->h, a[j], &absent);
//                     kh_key(g_pos->h, k) = (pos[j]&YAK_POS_MASK);
// 					// 输出kmer和pos
//         			printf("insert kmer: 0x%016llx, pos: %u\n", (unsigned long long)a[j], pos[j]);
//                     ++n_ins;
//                 }else if(((kh_key(g->h, k)&YAK_KEY_MASK) ^ r[j]) != 0){
//                     kh_key(g->h,k) |= YAK_REPEAT_MASK;
//                 }
//             }
//         } else {
//             k = yak_ht_get(g->h, (a[j] >> h->pre)<<YAK_COUNTER_BITS);
//             if(k != kh_end(g->h) && ((kh_key(g->h, k)&YAK_KEY_MASK) ^ r[j]) != 0){
//                 kh_key(g->h,k) |= YAK_REPEAT_MASK;
//             }
//         }
//     }
//     return n_ins;
// }

int yak_ch_insert_list_kmer_pos2(yak_ch_t *h, yak_ch_t *h_pos, int create_new, int n, const uint64_t *a, const uint16_t *r, const uint32_t *pos, long i)
{
	int j, mask = (1 << h->pre) - 1, n_ins = 0;
	yak_ch1_t *g = NULL;
	yak_ch1_t *g_pos = NULL;
	if (n == 0)
		return 0;
	g = &h->h[a[0] & mask];
	g_pos = &h_pos->h[a[0] & ((1 << h_pos->pre) - 1)];

	for (j = 0; j < n; ++j)
	{
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		// 只处理主表h中已存在的kmer
		if ((a[j] & mask) != (a[0] & mask))
			continue;
		k = yak_ht_get(g->h, x << YAK_COUNTER_BITS);
		if (k == kh_end(g->h))
			continue; // 主表不存在则跳过

		// 插入或更新h_pos
		int absent;
		khint_t k_pos = yak_ht_put(g_pos->h, a[j], &absent);
		kh_key(g_pos->h, k_pos) = (pos[j] & YAK_POS_MASK);

		std::string kmer_str = decode_kmer(a[j], h->k);
		// printf("insert h_pos: kmer: %s, pos: %u\n", kmer_str.c_str(), pos[j]);
		++n_ins;
		// 输出kmer和pos
		// printf("insert h_pos: 0x%016llx, pos: %u\n", (unsigned long long)a[j], pos[j]);
		++n_ins;
	}
	return n_ins;
}
int yak_ch_insert_list_kmer_pos(yak_ch_t *h, yak_ch_t *h_pos, int create_new, int n, const uint64_t *a, const uint16_t *r, const uint32_t *pos, long i)
{
	int j, mask = (1 << h->pre) - 1, n_ins = 0;
	int mask_pos = (1 << 30) - 1;
	yak_ch1_t *g;
	yak_ch1_t *g_pos;
	// recordset_ps_t *cur_recordset = &recordset[a[0]&mask];
	if (n == 0)
		return 0;
	g = &h->h[a[0] & mask];
	g_pos = &h_pos->h[a[0] & mask_pos];
	for (j = 0; j < n; ++j)
	{
		int ins = 1, absent;
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		if ((a[j] & mask) != (a[0] & mask))
			continue;
		if (create_new)
		{
			if (g->b)
				ins = (yak_bf_insert(g->b, (a[j] >> h->pre)) == h->n_hash);
			if (ins)
			{
				k = yak_ht_put(g->h, (a[j] >> h->pre) << YAK_COUNTER_BITS, &absent);
				if (absent)
				{
					kh_key(g->h, k) |= r[j];
					// --- 解析 h 中信息 ---
					uint64_t key_h = kh_key(g->h, k);
					uint64_t kmer_code = key_h >> YAK_COUNTER_BITS;
					uint64_t read_id = key_h & YAK_KEY_MASK;
					bool is_forward = key_h & YAK_FORWARD_MASK;
					bool is_repeat = key_h & YAK_REPEAT_MASK;

					std::string kmer_str = decode_kmer(a[j], h->k);
					printf("from h: kmer=%s, readID=%lu, forward=%d, repeat=%d\n",
						   kmer_str.c_str(), read_id, is_forward, is_repeat);

					k = yak_ht_put(g_pos->h, (a[j] >> 30) << 30, &absent);
					kh_key(g_pos->h, k) |= (pos[j] & YAK_POS_MASK);
					++n_ins;

					// --- 解析 h_pos 中信息 ---
					uint64_t key_pos = kh_key(g_pos->h, k);
					uint32_t position = key_pos & YAK_POS_MASK;
					printf("from h_pos: kmer_high30=0x%08llx, pos=%u\n",
						   (unsigned long long)(a[j] >> 30), position);
				}
				else if (((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
				{
					kh_key(g->h, k) |= YAK_REPEAT_MASK;
				}
			}
		}
		else
		{
			k = yak_ht_get(g->h, (a[j] >> h->pre) << YAK_COUNTER_BITS);
			if (k != kh_end(g->h) && ((kh_key(g->h, k) & YAK_KEY_MASK) ^ r[j]) != 0)
			{
				kh_key(g->h, k) |= YAK_REPEAT_MASK;
			}
		}
	}
	return n_ins;
}

// int yak_ch_insert_list_kmer_record_mapping(
//     yak_ch_t *h, int create_new, int n, const uint64_t *a, const bool *f, const int *pos, const uint16_t *r, long i, std::set<uint64_t>* deleted
// )
// {
//     int j, mask = (1<<h->pre) - 1, n_ins = 0;
//     yak_ch1_t *g;
//     if (n == 0) return 0;
//     g = &h->h[a[0]&mask];
//     for (j = 0; j < n; ++j) {
//         int ins = 1, absent;
//         uint64_t x = a[j] >> h->pre;
//         khint_t k;
//         if ((a[j]&mask) != (a[0]&mask)) continue;
//         if (create_new) {
//             if (g->b)
//                 ins = (yak_bf_insert(g->b, x) == h->n_hash);
//             if (ins) {
//                 k = yak_ht_put(g->h, x<<YAK_COUNTER_BITS, &absent);
//                 if (absent){
//                     ++n_ins;
//                     kh_key(g->h, k)|=r[j];
//                     if(f[j]){
//                         kh_key(g->h, k) |= YAK_FORWARD_MASK;
//                     }
//                     // 这里可以处理pos，比如打印、存储到另一个表、或编码到key的高位（需保证不会冲突）
//                     // 例如：你可以在此处调用一个自定义函数保存pos信息
//                     // save_kmer_pos(x, pos[j]);
//                 }else if(((kh_key(g->h, k)&YAK_KEY_MASK) ^ r[j]) != 0){
//                     kh_key(g->h,k) |= YAK_REPEAT_MASK;
//                 }
//             }
//         } else {
//             k = yak_ht_get(g->h, x<<YAK_COUNTER_BITS);
//             if(k != kh_end(g->h) && ((kh_key(g->h, k)&YAK_KEY_MASK) ^ r[j]) != 0){
//                 kh_key(g->h,k) |= YAK_REPEAT_MASK;
//             }
//         }
//     }
//     return n_ins;
// }
uint16_t yak_ch_get_pos2(const yak_ch_t *h, const yak_ch_t *h_pos, uint64_t x, uint64_t x2, uint32_t *pos)
{
	int mask = (1 << h->pre) - 1;
	int mask_pos = (1 << 30) - 1;
	yak_ht_t *g = h->h[x & mask].h;
	yak_ht_t *g_pos = h_pos->h[x & mask_pos].h;
	khint_t k, k2;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	if (k == kh_end(g) || (kh_key(g, k) & YAK_REPEAT_MASK) != 0)
	{
		return 65535;
	}
	uint16_t id = kh_key(g, k) & YAK_KEY_MASK;
	k = yak_ht_get(g_pos, (x >> 30) << 30);
	k2 = yak_ht_get(g_pos, (x2 >> 30) << 30);
	if (k != kh_end(g_pos))
	{
		(*pos) = (kh_key(g_pos, k) & YAK_POS_MASK);
		printf("use k1 pos\n");
	}

	else if (k2 != kh_end(g_pos))
	{
		(*pos) = (kh_key(g_pos, k2) & YAK_POS_MASK);
		printf("use k2 pos\n");
	}
	return id;
}

// void yak_ch_get_pos_test(const yak_ch_t *h, const yak_ch_t *h_pos, uint64_t x, uint32_t *pos)
// {
// 	khint_t k;
// 	int mask_pos = (1 << 30) - 1;
// 	yak_ht_t *g_pos = h_pos->h[x & mask_pos].h;
// 	k = yak_ht_get(g_pos, (x >> 30) << 30);
// 	uint64_t key_pos = kh_key(g_pos, k);
// 	uint32_t position = key_pos & YAK_POS_MASK;
// 	if (key_pos != kh_end(g_pos))
// 	{
// 		(*pos) = position;
// 		printf("position_test=%u\n", position);
// 	}
// }

void yak_ch_get_pos_test(yak_ch_t *h, yak_ch_t *h_pos, uint64_t x, uint32_t *pos)
{
	int mask_pos = (1 << 30) - 1;
	yak_ht_t *g_pos = h_pos->h[x & mask_pos].h;

	khint_t k = yak_ht_get(g_pos, (x >> 30) << 30); // 或者 (x >> h_pos->pre) << h_pos->pre
	if (k == kh_end(g_pos))
	{

		*pos = 65535;
		return;
	}

	uint64_t key_pos = kh_key(g_pos, k);
	uint32_t position = key_pos & YAK_POS_MASK;
	*pos = position;
	printf("position_test=%u\n", position);
}

uint64_t yak_ch_get_pos(const yak_ch_t *h, const yak_ch_t *h_pos, uint64_t x, uint32_t *pos)
{
	int mask = (1 << h->pre) - 1;
	int mask_pos = (1 << 30) - 1;
	yak_ht_t *g = h->h[x & mask].h;
	yak_ht_t *g_pos = h_pos->h[x & mask_pos].h;
	// std::cout<<"yak_ch_get_pos : x & mask_pos :" << (x & mask_pos) << std::endl;
	khint_t k, k_pos;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	if (k == kh_end(g) || (kh_key(g, k) & YAK_REPEAT_MASK) != 0)
	{
		return -1;
	}

	uint64_t id = kh_key(g, k) & YAK_KEY_MASK;
	// std::cout<<(x >> 30) << 30<<std::endl;
	k_pos = yak_ht_get(g_pos, (x >> 30) << 30);
	// std::cout<<"yak_ch_get_pos : (x >> 30) << 30 :" << ((x >> 30) << 30) << std::endl;
	if (k_pos != kh_end(g_pos))
	{
		(*pos) = (kh_key(g_pos, k_pos) & YAK_POS_MASK);
		// printf("Found position: %u\n", *pos); // 正常找到
		// std::cout<<"g_pos->count: "<<g_pos->count<<" g_pos->keys: "<<g_pos->keys<<std::endl;
	}
	else
	{

		if (k_pos == kh_end(g_pos))
		{
			// std::cout << k_pos << " " << kh_end(g_pos) << std::endl;
			// std::cout << "Position not found for x=" << x << std::endl;
			// std::cout << "g_pos->count: " << g_pos->count << " g_pos->keys: " << g_pos->keys << std::endl;

			(*pos) = -1;
			// uint32_t pos2 = kh_key(g_pos, k_pos) & YAK_KEY_MASK;
			// printf("pos2: %u\n", pos2);
		}
	}

	return id;
}

uint64_t yak_ch_get_pos_for_repeat(const yak_ch_t *h, const yak_ch_t *h_pos, uint64_t x, uint32_t *pos)
{
	int mask = (1 << h->pre) - 1;
	int mask_pos = (1 << 30) - 1;
	yak_ht_t *g = h->h[x & mask].h;
	yak_ht_t *g_pos = h_pos->h[x & mask_pos].h;
	// std::cout<<"yak_ch_get_pos : x & mask_pos :" << (x & mask_pos) << std::endl;
	khint_t k, k_pos;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	if (k == kh_end(g))
	{
		return -1;
	}
	uint64_t id = kh_key(g, k) & YAK_KEY_MASK;
	// std::cout<<(x >> 30) << 30<<std::endl;
	k_pos = yak_ht_get(g_pos, (x >> 30) << 30);
	// std::cout<<"yak_ch_get_pos : (x >> 30) << 30 :" << ((x >> 30) << 30) << std::endl;
	if (k_pos != kh_end(g_pos))
	{
		(*pos) = (kh_key(g_pos, k_pos) & YAK_POS_MASK);
		// printf("Found position: %u\n", *pos); // 正常找到
		//  std::cout<<"g_pos->count: "<<g_pos->count<<" g_pos->keys: "<<g_pos->keys<<std::endl;
	}
	else
	{

		if (k_pos == kh_end(g_pos))
		{
			std::cout << k_pos << " " << kh_end(g_pos) << std::endl;
			// std::cout << "Position not found for x=" << x << std::endl;
			// std::cout << "g_pos->count: " << g_pos->count << " g_pos->keys: " << g_pos->keys << std::endl;

			(*pos) = -1;
			// uint32_t pos2 = kh_key(g_pos, k_pos) & YAK_KEY_MASK;
			// printf("pos2: %u\n", pos2);
		}
	}

	return id;
}

int yak_ch_get(const yak_ch_t *h, uint64_t x)
{
	int mask = (1 << h->pre) - 1;
	yak_ht_t *g = h->h[x & mask].h;
	khint_t k;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	return k == kh_end(g) ? -1 : kh_key(g, k) & YAK_MAX_COUNT;
}

int yak_ch__non_repeat_get_k(const yak_ch_t *h, uint64_t x)
{
	int mask = (1 << h->pre) - 1;
	yak_ht_t *g = h->h[x & mask].h;
	khint_t k;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	return (k == kh_end(g)) ? -1 : kh_key(g, k);
}

// int yak_ch_get_full_info(const yak_ch_t *h, const yak_ch_t *h_pos, uint64_t x,
//                                        uint64_t *read_id, bool *is_forward, bool *is_repeat,
//                                        uint32_t *position)
// {
//     int mask = (1 << h->pre) - 1;
//     int mask_pos = (1 << 30) - 1;

//     // 查找 h 中的 k-mer key
//     yak_ht_t *g = h->h[x & mask].h;
//     khint_t k = yak_ht_get(g, (x >> h->pre) << YAK_COUNTER_BITS);
//     if (k == kh_end(g) || (kh_key(g, k) & YAK_REPEAT_MASK) != 0) {
//         return -1; // k-mer 不存在或者重复
//     }

//     uint64_t key_h = kh_key(g, k);
//     uint64_t kmer = key_h >> YAK_COUNTER_BITS;      // k-mer 本体
//     *read_id = key_h & YAK_KEY_MASK;               // read id
//     *is_forward = key_h & YAK_FORWARD_MASK;        // forward 标志
//     *is_repeat  = key_h & YAK_REPEAT_MASK;         // repeat 标志

//     // 查找 h_pos 中对应的位置
//     yak_ht_t *g_pos = h_pos->h[x & mask_pos].h;
//     uint64_t kmer_key_high = kmer >> 30 << 30; // 高位匹配 h_pos
//     khint_t k_pos = yak_ht_get(g_pos, kmer_key_high);
//     if (k_pos != kh_end(g_pos)) {
//         uint64_t key_pos = kh_key(g_pos, k_pos);
//         *position = key_pos & YAK_POS_MASK;
//     } else {
//         *position = 0; // 如果没找到，返回0
//     }

//     return 0; // 成功
// }
// 获取 kmer 对应的 read id、forward/repeat 标记，以及位置信息
// 返回 read id（如果不存在返回 65535），pos 通过指针返回
uint16_t yak_ch_get_full_info(const yak_ch_t *h, const yak_ch_t *h_pos, uint64_t x,
							  bool *is_forward, bool *is_repeat, uint32_t *position)
{
	int mask = (1 << h->pre) - 1;
	int mask_pos = (1 << 30) - 1;

	// 查找 kmer 本体
	yak_ht_t *g = h->h[x & mask].h;
	khint_t k = yak_ht_get(g, (x >> h->pre) << YAK_COUNTER_BITS);
	// if (k == kh_end(g) || (kh_key(g, k) & YAK_REPEAT_MASK)) {
	//     // kmer 不存在或者重复
	//     if (position) *position = 0;
	//     if (is_forward) *is_forward = false;
	//     if (is_repeat) *is_repeat = true;
	//     return 65535;
	// }

	uint64_t key_h = kh_key(g, k);
	uint16_t read_id = key_h & YAK_KEY_MASK;
	if (is_forward)
		*is_forward = key_h & YAK_FORWARD_MASK;
	if (is_repeat)
		*is_repeat = key_h & YAK_REPEAT_MASK;

	// 查找位置信息
	yak_ht_t *g2 = h_pos->h[x & mask_pos].h;
	yak_ch1_t *g_pos = &h_pos->h[x & mask_pos];
	khint_t k_pos = yak_ht_get(g_pos->h, (x >> 30) << 30);
	if (k_pos != kh_end(g_pos->h) && position)
	{
		*position = kh_key(g_pos->h, k_pos) & YAK_POS_MASK;
		// printf("Found position: %u\n", *position); // 正常找到
	}
	// else if (position)
	// {
	// 	*position = 0;
	// 	if (k_pos == kh_end(g_pos->h))
	// 	{
	// 		if ((k == kh_end(g) || (kh_key(g, k) & YAK_REPEAT_MASK) != 0) ? -1 : kh_key(g, k) != -1)
	// 			printf("Position not found in h_pos table,but kmer found in h \n");
	// 	}
	// 	else
	// 	{
	// 		printf("Position pointer is valid but key exists, setting position=0\n");
	// 	}
	// }

	return (k == kh_end(g) || (kh_key(g, k) & YAK_REPEAT_MASK) != 0) ? -1 : kh_key(g, k);
}

int yak_ch_get_k(const yak_ch_t *h, uint64_t x)
{
	int mask = (1 << h->pre) - 1;
	yak_ht_t *g = h->h[x & mask].h;
	khint_t k;
	// std::cout<<"g : x & mask :" << (x & mask) << std::endl;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	// if ((kh_key(g, k) & YAK_REPEAT_MASK) != 0)
	// {
	// 	std::cout << "kmer is repeat" << std::endl;
	// }
	return (k == kh_end(g) || (kh_key(g, k) & YAK_REPEAT_MASK) != 0) ? -1 : kh_key(g, k);
}

int yak_ch_get_k_for_repeat(const yak_ch_t *h, uint64_t x)
{
	int mask = (1 << h->pre) - 1;
	yak_ht_t *g = h->h[x & mask].h;
	khint_t k;
	// std::cout<<"g : x & mask :" << (x & mask) << std::endl;
	k = yak_ht_get(g, x >> h->pre << YAK_COUNTER_BITS);
	// if ((kh_key(g, k) & YAK_REPEAT_MASK) != 0)
	// {
	// 	std::cout << "kmer is repeat" << std::endl;
	// }
	return (k == kh_end(g)) ? -1 : kh_key(g, k);
}
/*************************
 * Clear all counts to 0 *
 *************************/

static void worker_clear(void *data, long i, int tid) // callback for kt_for()
{
	yak_ch_t *h = (yak_ch_t *)data;
	yak_ht_t *g = h->h[i].h;
	khint_t k;
	uint64_t mask = ~1ULL >> YAK_COUNTER_BITS << YAK_COUNTER_BITS;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			kh_key(g, k) &= mask;
}

void yak_ch_clear(yak_ch_t *h, int n_thread)
{
	kt_for(n_thread, worker_clear, h, 1 << h->pre);
}

/*************
 * Histogram *
 *************/

typedef struct
{
	uint64_t c[YAK_N_COUNTS];
} buf_cnt_t;

typedef struct
{
	const yak_ch_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

static void worker_hist(void *data, long i, int tid) // callback for kt_for()
{
	hist_aux_t *a = (hist_aux_t *)data;
	uint64_t *cnt = a->cnt[tid].c;
	yak_ht_t *g = a->h->h[i].h;
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			++cnt[kh_key(g, k) & YAK_MAX_COUNT];
}

void yak_ch_hist(const yak_ch_t *h, int64_t cnt[YAK_N_COUNTS], int n_thread)
{
	hist_aux_t a;
	int i, j;
	a.h = h;
	memset(cnt, 0, YAK_N_COUNTS * sizeof(uint64_t));
	CALLOC(a.cnt, n_thread);
	kt_for(n_thread, worker_hist, &a, 1 << h->pre);
	for (i = 0; i < YAK_N_COUNTS; ++i)
		cnt[i] = 0;
	for (j = 0; j < n_thread; ++j)
		for (i = 0; i < YAK_N_COUNTS; ++i)
			cnt[i] += a.cnt[j].c[i];
	free(a.cnt);
}

/**********
 * Shrink *
 **********/

typedef struct
{
	int min, max;
	yak_ch_t *h;
} shrink_aux_t;

static void worker_shrink(void *data, long i, int tid) // callback for kt_for()
{
	shrink_aux_t *a = (shrink_aux_t *)data;
	yak_ch_t *h = a->h;
	yak_ht_t *g = h->h[i].h, *f;
	khint_t k;
	f = yak_ht_init();
	yak_ht_resize(f, kh_size(g));
	for (k = 0; k < kh_end(g); ++k)
	{
		if (kh_exist(g, k))
		{
			int absent, c = kh_key(g, k) & YAK_MAX_COUNT;
			if (c >= a->min && c <= a->max)
				yak_ht_put(f, kh_key(g, k), &absent);
		}
	}
	yak_ht_destroy(g);
	h->h[i].h = f;
}

void yak_ch_shrink(yak_ch_t *h, int min, int max, int n_thread)
{
	int i;
	shrink_aux_t a;
	a.h = h, a.min = min, a.max = max;
	kt_for(n_thread, worker_shrink, &a, 1 << h->pre);
	for (i = 0, h->tot = 0; i < 1 << h->pre; ++i)
		h->tot += kh_size(h->h[i].h);
}

/*******
 * I/O *
 *******/

int yak_ch_dump(const yak_ch_t *ch, const char *fn)
{
	FILE *fp;
	uint32_t t[3];
	int i;
	if ((fp = strcmp(fn, "-") ? fopen(fn, "wb") : stdout) == 0)
		return -1;
	fwrite(YAK_MAGIC, 1, 4, fp);
	t[0] = ch->k, t[1] = ch->pre, t[2] = YAK_COUNTER_BITS;
	fwrite(t, 4, 3, fp);
	for (i = 0; i < 1 << ch->pre; ++i)
	{
		yak_ht_t *h = ch->h[i].h;
		khint_t k;
		t[0] = kh_capacity(h), t[1] = kh_size(h);
		fwrite(t, 4, 2, fp);
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k))
				fwrite(&kh_key(h, k), 8, 1, fp);
	}
	fprintf(stderr, "[M::%s] dumpped the hash table to file '%s'.\n", __func__, fn);
	fclose(fp);
	return 0;
}

yak_ch_t *yak_ch_restore_core(yak_ch_t *ch0, const char *fn, int mode, ...)
{
	va_list ap;
	FILE *fp;
	uint32_t t[3];
	char magic[4];
	int i, j, absent, min_cnt = 0, mid_cnt = 0, mode_err = 0;
	uint64_t mask = (1ULL << YAK_COUNTER_BITS) - 1, n_ins = 0, n_new = 0;
	yak_ch_t *ch;

	va_start(ap, mode);
	if (mode == YAK_LOAD_ALL)
	{ // do nothing
	}
	else if (mode == YAK_LOAD_TRIOBIN1 || mode == YAK_LOAD_TRIOBIN2)
	{
		assert(YAK_COUNTER_BITS >= 4);
		min_cnt = va_arg(ap, int);
		mid_cnt = va_arg(ap, int);
		if (ch0 == 0 && mode == YAK_LOAD_TRIOBIN2)
			mode_err = 1;
	}
	else
		mode_err = 1;
	va_end(ap);
	if (mode_err)
		return 0;

	if ((fp = fopen(fn, "rb")) == 0)
		return 0;
	if (fread(magic, 1, 4, fp) != 4)
		return 0;
	if (strncmp(magic, YAK_MAGIC, 4) != 0)
	{
		fprintf(stderr, "ERROR: wrong file magic.\n");
		fclose(fp);
		return 0;
	}
	fread(t, 4, 3, fp);
	if (t[2] != YAK_COUNTER_BITS)
	{
		fprintf(stderr, "ERROR: saved counter bits: %d; compile-time counter bits: %d\n", t[2], YAK_COUNTER_BITS);
		fclose(fp);
		return 0;
	}

	ch = ch0 == 0 ? yak_ch_init(t[0], t[1], 0, 0) : ch0;
	assert((int)t[0] == ch->k && (int)t[1] == ch->pre);
	for (i = 0; i < 1 << ch->pre; ++i)
	{
		yak_ht_t *h = ch->h[i].h;
		fread(t, 4, 2, fp);
		if (ch0 == 0)
			yak_ht_resize(h, t[0]);
		for (j = 0; j < t[1]; ++j)
		{
			uint64_t key;
			fread(&key, 8, 1, fp);
			if (mode == YAK_LOAD_ALL)
			{
				++n_ins;
				yak_ht_put(h, key, &absent);
				if (absent)
					++n_new;
			}
			else if (mode == YAK_LOAD_TRIOBIN1 || mode == YAK_LOAD_TRIOBIN2)
			{
				int cnt = key & mask, x, shift = mode == YAK_LOAD_TRIOBIN1 ? 0 : 2;
				if (cnt >= mid_cnt)
					x = 2 << shift;
				else if (cnt >= min_cnt)
					x = 1 << shift;
				else
					x = -1;
				if (x >= 0)
				{
					khint_t k;
					key = (key & ~mask) | x;
					++n_ins;
					k = yak_ht_put(h, key, &absent);
					if (absent)
						++n_new;
					else
						kh_key(h, k) = kh_key(h, k) | x;
				}
			}
		}
	}
	fclose(fp);
	fprintf(stderr, "[M::%s] inserted %ld k-mers, of which %ld are new\n", __func__, (long)n_ins, (long)n_new);
	return ch;
}

yak_ch_t *yak_ch_restore(const char *fn)
{
	return yak_ch_restore_core(0, fn, YAK_LOAD_ALL);
}
