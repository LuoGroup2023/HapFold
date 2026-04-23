extern "C"
{
#include "paf.h"
}
#include <string.h>
#include <iostream>
#include "graph.h"
#include "sys.cpp"
#include <fstream>
#include "bubble_chain.h"
#include "resolve_repeat_haplotype.h"
#include "hic_mapping.h"
#define HapFold_VERSION "0.1"

typedef struct
{
	int num_reads;
	char **read_file_names;
	char *output_file_name;
	int thread_num;
} ps_opt_t;

ps_opt_t asm_opt;

void init_opt(ps_opt_t *asm_opt)
{
	memset(asm_opt, 0, sizeof(ps_opt_t));
	asm_opt->num_reads = 0;
	asm_opt->read_file_names = NULL;
	asm_opt->thread_num = 1;
}

void destory_opt(ps_opt_t *asm_opt)
{
	if (asm_opt->read_file_names != NULL)
	{
		free(asm_opt->read_file_names);
	}
}




static ko_longopt_t long_options[] = {
	{"version", ko_no_argument, 300},
	{"write-paf", ko_no_argument, 302},
	{0, 0, 0}};

void Print_H(ps_opt_t *asm_opt)
{
	fprintf(stderr, "Usage: HapFold [options] <input> <output> <...>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  Assembly:\n");
	fprintf(stderr, "    -o FILE       prefix of output files [%s]\n", asm_opt->output_file_name);
	fprintf(stderr, "    -t INT        number of threads [%d]\n", asm_opt->thread_num);
	fprintf(stderr, "    --version     show version number\n");
	fprintf(stderr, "    -h            show help information\n");

	fprintf(stderr, "Example: ./HapFold bubble_chain myfile.gfa -o []\n");
	fprintf(stderr, "Example: ./HapFold resolve_repeats myfile.gfa [].paf.gz -o []\n");
	fprintf(stderr, "Example: ./HapFold asm_component \n");
	fprintf(stderr, "Example: ./HapFold filter_paf -t[n] myfile.gfa [].paf.gz -o [].paf \n");
	//fprintf(stderr, "See `man ./HapFold.1' for detailed description of these command-line options.\n");
}

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g')
		x *= 1e9;
	else if (*p == 'M' || *p == 'm')
		x *= 1e6;
	else if (*p == 'K' || *p == 'k')
		x *= 1e3;
	return (int64_t)(x + .499);
}

std::vector<NamedBubbleContig> read_named_bubble_contigs(const std::string &filename)
{
	std::ifstream infile(filename);
	std::vector<NamedBubbleContig> result;

	if (!infile)
	{
		std::cerr << "[ERROR] Cannot open file: " << filename << "\n";
		return result;
	}

	std::string line;
	std::vector<std::string> buffer;

	while (std::getline(infile, line))
	{
		if (line.empty())
			continue;

		std::vector<std::string> names;
		std::stringstream ss(line);
		std::string token;
		while (std::getline(ss, token, ','))
		{
			if (!token.empty())
				names.push_back(token);
		}

		buffer.push_back(line);

		if (buffer.size() == 2)
		{
			NamedBubbleContig contig;
			std::stringstream ss1(buffer[0]), ss2(buffer[1]);
			std::string name;

			while (std::getline(ss1, name, ','))
			{
				if (!name.empty())
					contig.hap1.push_back(name);
			}
			while (std::getline(ss2, name, ','))
			{
				if (!name.empty())
					contig.hap2.push_back(name);
			}

			result.push_back(std::move(contig));
			buffer.clear();
		}
	}

	if (!buffer.empty())
	{
		std::cerr << "[WARNING] File " << filename << " has odd number of lines, last bubble is incomplete.\n";
	}

	return result;
}
int main_resolve_haplotypes(int argc, char *argv[])
{
    ketopt_t o = KETOPT_INIT;
    int c;

    GlobalParams g_params;

    // 1. 增加控制 12M 和 300K 的长参数选项 (代号 303, 304)
    static ko_longopt_t longopts[] = {
        {"hic_scaffold_threshold_ratio", ko_required_argument, 301},
        {"debug", ko_no_argument, 302}, 
        {"chain_len_thresh", ko_required_argument, 303},     // 对应 > 12M 参与迭代的阈值
        {"scaffold_len_thresh", ko_required_argument, 304},  // 对应 > 300K 直接输出的阈值
        {0, 0, 0} 
    };

    while ((c = ketopt(&o, argc, argv, 1, "t:e:i:f:1:2:u:c:n:pd", longopts)) >= 0)
    {
        if (c == 't')
            g_params.n_threads = atoi(o.arg);
        else if (c == 'e')
            g_params.enzymes_unsplit = string(o.arg);
        else if (c == 'i')
            g_params.check_identity = (strcmp(o.arg, "true") == 0);
        else if (c == 'f')
            g_params.identityFile = string(o.arg);
        else if (c == '1')
            g_params.hap1_gfa = string(o.arg);
        else if (c == '2')
            g_params.hap2_gfa = string(o.arg);
        else if (c == 'u')
            g_params.utg_ctg_file = string(o.arg);
        else if (c == 'c')
            g_params.contig_hap_file = string(o.arg);
        else if (c == 'n')
            g_params.n_chrs = atoi(o.arg);
        else if (c == 'p') 
            g_params.is_plant = true;
        else if (c == 'd' || c == 302) 
            g_params.debug_mode = true;
        else if (c == 301) 
            g_params.hic_scaffold_threshold_ratio = atof(o.arg);
        else if (c == 303) // 捕获 12M 阈值
            g_params.chain_len_threshold = atoi(o.arg);
        else if (c == 304) // 捕获 300K 阈值
            g_params.scaffold_len_threshold = atoi(o.arg);
    }

    if (argc - o.ind < 3)
    {
        fprintf(stderr, "\nUsage: GraPhaser resolve_haplotypes [options] <hic_mapping.out> <assembly.gfa> <output_dir> -1 *.hap1.p_ctg.gfa -2 *.hap2.p_ctg.gfa -u utg_ctg_file \n\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -t INT      Number of threads [%d]\n", g_params.n_threads);
        fprintf(stderr, "  -n INT      Expected number of chromosomes (e.g., 78 for chicken) [%d]\n", g_params.n_chrs);
        fprintf(stderr, "  -e STR      Restriction enzymes separated by comma (e.g., GATC,GANTC) [%s]\n", g_params.enzymes_unsplit.c_str());
        fprintf(stderr, "  -c FILE     Path to contig_hap_nodes.txt (debug for Hi-C phasing)\n");
        fprintf(stderr, "  -u FILE     Path to utg_to_ctg relationship file [optional]\n");
        fprintf(stderr, "  -1 FILE     Path to haplotype 1 GFA file (*.hap1.p_ctg.gfa)\n");
        fprintf(stderr, "  -2 FILE     Path to haplotype 2 GFA file (*.hap2.p_ctg.gfa)\n");
        fprintf(stderr, "  -i BOOL     Enable identity check on contigs (true/false) [%s]\n", (g_params.check_identity ? "true" : "false"));
        fprintf(stderr, "  -f FILE     Precomputed identity file path; if omitted, check will run automatically [%s]\n", g_params.identityFile.c_str());
        fprintf(stderr, "  -p          Enable plant mode (uses alternative phasing algorithms) [optional]\n"); 
        fprintf(stderr, "  -d, --debug Enable debug mode to run test code functions [optional]\n"); 
        fprintf(stderr, "  --hic_scaffold_threshold_ratio FLOAT  Threshold ratio for Hi-C scaffolding [%.2f]\n", g_params.hic_scaffold_threshold_ratio);
        fprintf(stderr, "  --chain_len_thresh INT                Length threshold to join contig_chain for iterative merging [%d]\n", g_params.chain_len_threshold);
        fprintf(stderr, "  --scaffold_len_thresh INT             Length threshold to directly output to scaffold.fa [%d]\n", g_params.scaffold_len_threshold);
        fprintf(stderr, "\n");
        return 1;
    }

    // vector<string> enzymes;
    // if (g_params.enzymes_unsplit.size() > 1)
    // {
    //     stringstream s_stream(g_params.enzymes_unsplit);
    //     while (s_stream.good())
    //     {
    //         string substr;
    //         getline(s_stream, substr, ',');
    //         substr.erase(remove(substr.begin(), substr.end(), '^'), substr.end());
    //         enzymes.push_back(substr);
    //     }
    // }
    std::vector<NamedBubbleContig> named_bubble_contigs;

    if (!g_params.contig_hap_file.empty())
    {
        named_bubble_contigs = read_named_bubble_contigs(g_params.contig_hap_file);
        std::cerr << "[INFO] Loaded " << named_bubble_contigs.size() << " named bubble contigs from " << g_params.contig_hap_file << "\n";
    }

    char *connectionFile = argv[o.ind];
    char *gfa_filename = argv[o.ind + 1];
    char *output_directory = argv[o.ind + 2];
    printf("start main\n");
    asg_t *graph = gfa_read(gfa_filename);
    map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_graph = nullptr;
    
    uint32_t **connections_foward;
    CALLOC(connections_foward, graph->n_seq);
    for (int i = 0; i < graph->n_seq; i++)
    {
        CALLOC(connections_foward[i], graph->n_seq);
        memset(connections_foward[i], 0, sizeof(*connections_foward[i]));
    }
    uint32_t **connections_backward;
    CALLOC(connections_backward, graph->n_seq);
    for (int i = 0; i < graph->n_seq; i++)
    {
        CALLOC(connections_backward[i], graph->n_seq);
        memset(connections_backward[i], 0, sizeof(*connections_backward[i]));
    }
    ifstream infile(connectionFile);
    uint32_t i, j, count_forward, count_backward;
    while (infile >> i >> j >> count_backward >> count_forward)
    {
        connections_backward[i][j] = count_backward;
        connections_backward[j][i] = count_backward;
        connections_foward[i][j] = count_forward;
        connections_foward[j][i] = count_forward;
    }
    
    std::string utg_gfa = std::string(gfa_filename);

    if (g_params.is_plant)
    {
        printf("[INFO] Plant mode enabled. Using alternative phasing functions.\n");
        bubble_chain_graph = phasing_plant_version(graph, string(output_directory), connections_foward, connections_backward);
    }
    else
    {
        printf("[INFO] Default mode enabled. Using standard phasing functions.\n");
        bubble_chain_graph = phasing_10_7(graph, string(output_directory), connections_foward, connections_backward);

        // 2. 将之前散落的各种参数整合进 g_params 统一传递
        if (g_params.debug_mode) {
            printf("[INFO] Debug mode enabled. Executing get_haplotype_path_test_code...\n");
            get_haplotype_path_test_code(connections_foward, connections_backward, graph, bubble_chain_graph,
                                         output_directory, named_bubble_contigs, gfa_filename, g_params);
        } else {
            printf("[INFO] Executing standard model get_haplotype_path_now...\n");
            get_haplotype_path_now(connections_foward, connections_backward, graph, bubble_chain_graph,
                                   output_directory, named_bubble_contigs, gfa_filename, g_params);
        }
    }
    return 0;
}




int main_hic_mapping(int argc, char *argv[])
{
	return main_poreC_map_test(argc, argv);
}

int main_count(int argc, char *argv[])
{
	// yak_ch_t *h;
	// int c;
	// char *fn_out = 0;
	// yak_copt_t opt;
	// ketopt_t o = KETOPT_INIT;
	// yak_copt_init(&opt);
	// while ((c = ketopt(&o, argc, argv, 1, "k:p:K:t:b:H:o:", 0)) >= 0)
	// {
	// 	if (c == 'k')
	// 		opt.k = atoi(o.arg);
	// 	else if (c == 'p')
	// 		opt.pre = atoi(o.arg);
	// 	else if (c == 'K')
	// 		opt.chunk_size = atoi(o.arg);
	// 	else if (c == 't')
	// 		opt.n_thread = atoi(o.arg);
	// 	else if (c == 'b')
	// 		opt.bf_shift = atoi(o.arg);
	// 	else if (c == 'H')
	// 		opt.bf_n_hash = mm_parse_num(o.arg);
	// 	else if (c == 'o')
	// 		fn_out = o.arg;
	// }
	// if (argc - o.ind < 1)
	// {
	// 	fprintf(stderr, "Usage: HapFold count [options] <in.fa> [in.fa]\n");
	// 	fprintf(stderr, "Options:\n");
	// 	fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
	// 	fprintf(stderr, "  -p INT     prefix length [%d]\n", opt.pre);
	// 	fprintf(stderr, "  -b INT     set Bloom filter size to 2**INT bits; 0 to disable [%d]\n", opt.bf_shift);
	// 	fprintf(stderr, "  -H INT     use INT hash functions for Bloom filter [%d]\n", opt.bf_n_hash);
	// 	fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
	// 	fprintf(stderr, "  -o FILE    dump the count hash table to FILE []\n");
	// 	fprintf(stderr, "  -K INT     chunk size [100m]\n");
	// 	fprintf(stderr, "Note: -b37 is recommended for human reads\n");
	// 	return 1;
	// }
	// if (opt.pre < YAK_COUNTER_BITS)
	// {
	// 	fprintf(stderr, "ERROR: -p should be at least %d\n", YAK_COUNTER_BITS);
	// 	return 1;
	// }
	// if (opt.k >= 64)
	// {
	// 	fprintf(stderr, "ERROR: -k must be smaller than 64\n");
	// 	return 1;
	// }
	// else if (opt.k >= 32)
	// {
	// 	fprintf(stderr, "WARNING: counts are inexact if -k is greater than 31\n");
	// }
	// h = yak_count(argv[o.ind], &opt, 0);
	// if (opt.bf_shift > 0)
	// {
	// 	yak_ch_destroy_bf(h);
	// 	yak_ch_clear(h, opt.n_thread);
	// 	h = yak_count(argc - o.ind >= 2 ? argv[o.ind + 1] : argv[o.ind], &opt, h);
	// 	yak_ch_shrink(h, 2, YAK_MAX_COUNT, opt.n_thread);
	// 	fprintf(stderr, "[M::%s] %ld distinct k-mers after shrinking\n", __func__, (long)h->tot);
	// }

	// if (fn_out)
	// 	yak_ch_dump(h, fn_out);
	// yak_ch_destroy(h);
	return 0;
}

int main(int argc, char *argv[])
{
	extern double yak_realtime(void);
	extern double yak_cputime(void);
	extern void yak_reset_realtime(void);
	double t_start;
	int ret = 0, i;

	if (argc == 1)
	{
		fprintf(stderr, "Usage: HapFold <command> <arguments> <inputs>\n");
		fprintf(stderr, "Commands:\n");
		fprintf(stderr, "  resolve_haplotypes    use hic data to resolve haplotypes\n");
		fprintf(stderr, "  hic_mapping           map hic data to sequences in the graph\n");
		fprintf(stderr, "  count                 count k-mers\n");
		fprintf(stderr, "  version               print version number\n");
		return 1;
	}
	yak_reset_realtime();
	t_start = yak_realtime();
	if (strcmp(argv[1], "resolve_haplotypes") == 0)
		ret = main_resolve_haplotypes(argc - 1, argv + 1);
	else if (strcmp(argv[1], "hic_mapping") == 0)
		ret = main_hic_mapping(argc - 1, argv + 1);
	else if (strcmp(argv[1], "count") == 0)
		ret = main_count(argc - 1, argv + 1);
	else if (strcmp(argv[1], "version") == 0)
	{
		printf("HapFold: %s\n", HapFold_VERSION);
		return 0;
	}
	else
	{
		fprintf(stderr, "[E::%s] unknown command\n", __func__);
		return 1;
	}
	if (ret == 0)
	{
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, HapFold_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, yak_realtime() - t_start, yak_cputime());
	}
	return ret;
}
