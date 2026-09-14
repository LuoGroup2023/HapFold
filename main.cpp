// /*
//  * HapFold - A graph-based haplotype reconstruction framework
//  * Copyright (C) 2024 Yichen Li 
//  *
//  * This program is free software: you can redistribute it and/or modify
//  * it under the terms of the GNU General Public License as published by
//  * the Free Software Foundation, either version 3 of the License, or
//  * (at your option) any later version.
//  * ...
//  */


// extern "C"
// {
// #include "paf.h"
// }
// #include <string.h>
// #include <iostream>
// #include "graph.h"
// // #include "sys.cpp"
// #include <fstream>
// #include "graph_refining.h"
// #include "phasing2scaffolding.h"
// #include "mapping.h"
// #define HapFold_VERSION "0.1"

// typedef struct
// {
// 	int num_reads;
// 	char **read_file_names;
// 	char *output_file_name;
// 	int thread_num;
// } ps_opt_t;

// ps_opt_t asm_opt;

// void init_opt(ps_opt_t *asm_opt)
// {
// 	memset(asm_opt, 0, sizeof(ps_opt_t));
// 	asm_opt->num_reads = 0;
// 	asm_opt->read_file_names = NULL;
// 	asm_opt->thread_num = 1;
// }

// void destory_opt(ps_opt_t *asm_opt)
// {
// 	if (asm_opt->read_file_names != NULL)
// 	{
// 		free(asm_opt->read_file_names);
// 	}
// }




// static ko_longopt_t long_options[] = {
// 	{"version", ko_no_argument, 300},
// 	{"write-paf", ko_no_argument, 302},
// 	{0, 0, 0}};

// void Print_H(ps_opt_t *asm_opt)
// {
//     fprintf(stderr, "\nUsage: HapFold <command> [options]\n\n");
    
//     fprintf(stderr, "Commands:\n");
//     fprintf(stderr, "    mapping       Map Hi-C/Pore-C reads to the unitig sequences\n");
//     fprintf(stderr, "    scaffolding   Refine graph, phase haplotypes, and build scaffolds\n\n");
    
//     fprintf(stderr, "Global Options:\n");
//     fprintf(stderr, "    -o FILE       prefix of output files/directory [%s]\n", asm_opt->output_file_name);
//     fprintf(stderr, "    -t INT        number of threads [%d]\n", asm_opt->thread_num);
//     fprintf(stderr, "    --version     show version number\n");
//     fprintf(stderr, "    -h            show help information\n\n");

//     fprintf(stderr, "Examples:\n");
//     fprintf(stderr, "  Step 1. Hi-C/Pore-C Mapping:\n");
//     fprintf(stderr, "    ./HapFold mapping -t 32 -1 hic_1.fq.gz -2 hic_2.fq.gz -o mapping.txt utg.fa\n\n");
    
//     fprintf(stderr, "  Step 2. Scaffolding (Graph refining & Phasing):\n");
//     fprintf(stderr, "    ./HapFold scaffolding mapping.txt assembly.gfa out_dir -1 hap1.p_ctg.gfa -2 hap2.p_ctg.gfa -u utg_ctg.csv\n\n");
// }

// static inline int64_t mm_parse_num(const char *str)
// {
// 	double x;
// 	char *p;
// 	x = strtod(str, &p);
// 	if (*p == 'G' || *p == 'g')
// 		x *= 1e9;
// 	else if (*p == 'M' || *p == 'm')
// 		x *= 1e6;
// 	else if (*p == 'K' || *p == 'k')
// 		x *= 1e3;
// 	return (int64_t)(x + .499);
// }

// std::vector<NamedBubbleContig> read_named_bubble_contigs(const std::string &filename)
// {
// 	std::ifstream infile(filename);
// 	std::vector<NamedBubbleContig> result;

// 	if (!infile)
// 	{
// 		std::cerr << "[ERROR] Cannot open file: " << filename << "\n";
// 		return result;
// 	}

// 	std::string line;
// 	std::vector<std::string> buffer;

// 	while (std::getline(infile, line))
// 	{
// 		if (line.empty())
// 			continue;

// 		std::vector<std::string> names;
// 		std::stringstream ss(line);
// 		std::string token;
// 		while (std::getline(ss, token, ','))
// 		{
// 			if (!token.empty())
// 				names.push_back(token);
// 		}

// 		buffer.push_back(line);

// 		if (buffer.size() == 2)
// 		{
// 			NamedBubbleContig contig;
// 			std::stringstream ss1(buffer[0]), ss2(buffer[1]);
// 			std::string name;

// 			while (std::getline(ss1, name, ','))
// 			{
// 				if (!name.empty())
// 					contig.hap1.push_back(name);
// 			}
// 			while (std::getline(ss2, name, ','))
// 			{
// 				if (!name.empty())
// 					contig.hap2.push_back(name);
// 			}

// 			result.push_back(std::move(contig));
// 			buffer.clear();
// 		}
// 	}

// 	if (!buffer.empty())
// 	{
// 		std::cerr << "[WARNING] File " << filename << " has odd number of lines, last bubble is incomplete.\n";
// 	}

// 	return result;
// }
// int main_phasing_scaffolding(int argc, char *argv[])
// {
//     ketopt_t o = KETOPT_INIT;
//     int c;

//     GlobalParams g_params;
//     g_params.n_chrs = -1;
//     static ko_longopt_t longopts[] = {
//         {"hic_scaffold_threshold_ratio", ko_required_argument, 301},
//         {"debug", ko_no_argument, 302}, 
//         {"chain_len_thresh", ko_required_argument, 303},     // 对应 > 12M 参与迭代的阈值
//         {"scaffold_len_thresh", ko_required_argument, 304},  // 对应 > 300K 直接输出的阈值
//         {0, 0, 0} 
//     };

//     while ((c = ketopt(&o, argc, argv, 1, "t:e:i:f:1:2:u:c:n:pd", longopts)) >= 0)
//     {
//         if (c == 't')
//             g_params.n_threads = atoi(o.arg);
//         else if (c == 'e')
//             g_params.enzymes_unsplit = string(o.arg);
//         else if (c == 'i')
//             g_params.check_identity = (strcmp(o.arg, "true") == 0);
//         else if (c == 'f')
//             g_params.identityFile = string(o.arg);
//         else if (c == '1')
//             g_params.hap1_gfa = string(o.arg);
//         else if (c == '2')
//             g_params.hap2_gfa = string(o.arg);
//         else if (c == 'u')
//             g_params.utg_ctg_file = string(o.arg);
//         else if (c == 'c')
//             g_params.contig_hap_file = string(o.arg);
//         else if (c == 'n')
//             g_params.n_chrs = atoi(o.arg);
//         // else if (c == 'p') 
//         //     g_params.is_plant = true;
//         else if (c == 'd' || c == 302) 
//             g_params.debug_mode = true;
//         else if (c == 301) 
//             g_params.hic_scaffold_threshold_ratio = atof(o.arg);
//         else if (c == 303) // 捕获 12M 阈值
//             g_params.chain_len_threshold = atoi(o.arg);
//         else if (c == 304) // 捕获 300K 阈值
//             g_params.scaffold_len_threshold = atoi(o.arg);
//     }

//     if (argc - o.ind < 3)
//     {
//         fprintf(stderr, "\nUsage: HapFold scaffolding [options] <mapping.txt> <assembly.gfa> <output_dir> -1 *.hap1.p_ctg.gfa -2 *.hap2.p_ctg.gfa -n chr_number\n\n");
//         fprintf(stderr, "Options:\n");
//         fprintf(stderr, "  -t INT      Number of threads [%d]\n", g_params.n_threads);
//         fprintf(stderr, "  -n INT      Expected number of chromosomes (e.g., 78 for chicken) [%d]\n", g_params.n_chrs);
//         fprintf(stderr, "  -e STR      Restriction enzymes separated by comma (e.g., GATC,GANTC) [%s]\n", g_params.enzymes_unsplit.c_str());
//         fprintf(stderr, "  -c FILE     Path to contig_hap_nodes.txt (debug for Hi-C phasing)\n");
//         fprintf(stderr, "  -u FILE     Output path/name for the UTG-CTG mapping file [default: <output_dir>/utg_ctg_mappings.csv]\n");
//         fprintf(stderr, "  -1 FILE     Path to haplotype 1 GFA file (*.hap1.p_ctg.gfa)\n");
//         fprintf(stderr, "  -2 FILE     Path to haplotype 2 GFA file (*.hap2.p_ctg.gfa)\n");
//         fprintf(stderr, "  -i BOOL     Enable identity check on contigs (true/false) [%s]\n", (g_params.check_identity ? "true" : "false"));
//         fprintf(stderr, "  -f FILE     Precomputed identity file path; if omitted, check will run automatically [%s]\n", g_params.identityFile.c_str());
//         // fprintf(stderr, "  -p          Enable plant mode (uses alternative phasing algorithms) [optional]\n"); 
//         fprintf(stderr, "  -d, --debug Enable debug mode to run test code functions [optional]\n"); 
//         fprintf(stderr, "  --hic_scaffold_threshold_ratio FLOAT  Threshold ratio for Hi-C scaffolding [%.2f]\n", g_params.hic_scaffold_threshold_ratio);
//         fprintf(stderr, "  --chain_len_thresh INT                Length threshold to join contig_chain for iterative merging [%d]\n", g_params.chain_len_threshold);
//         fprintf(stderr, "  --scaffold_len_thresh INT             Length threshold to directly output to scaffold.fa [%d]\n", g_params.scaffold_len_threshold);
//         fprintf(stderr, "\n");
//         return 1;
//     }

//     // vector<string> enzymes;
//     // if (g_params.enzymes_unsplit.size() > 1)
//     // {
//     //     stringstream s_stream(g_params.enzymes_unsplit);
//     //     while (s_stream.good())
//     //     {
//     //         string substr;
//     //         getline(s_stream, substr, ',');
//     //         substr.erase(remove(substr.begin(), substr.end(), '^'), substr.end());
//     //         enzymes.push_back(substr);
//     //     }
//     // }
//     std::vector<NamedBubbleContig> named_bubble_contigs;

//     if (!g_params.contig_hap_file.empty())
//     {
//         named_bubble_contigs = read_named_bubble_contigs(g_params.contig_hap_file);
//         std::cerr << "[INFO] Loaded " << named_bubble_contigs.size() << " named bubble contigs from " << g_params.contig_hap_file << "\n";
//     }

//     char *connectionFile = argv[o.ind];
//     char *gfa_filename = argv[o.ind + 1];
//     char *output_directory = argv[o.ind + 2];

//     if (g_params.utg_ctg_file.empty())
//     {
//         std::string out_dir = std::string(output_directory);
//         if (!out_dir.empty() && out_dir.back() != '/')
//             out_dir += "/";

//         g_params.utg_ctg_file = out_dir + "utg_ctg_mappings.csv";

//         fprintf(stderr, "[INFO] No -u provided. Using default UTG-CTG mapping file: %s\n",
//                 g_params.utg_ctg_file.c_str());
//     }
//     else
//     {
//         fprintf(stderr, "[INFO] Using user-specified UTG-CTG mapping file: %s\n",
//                 g_params.utg_ctg_file.c_str());
//     }


//     if (g_params.hap1_gfa.empty() || g_params.hap2_gfa.empty() || g_params.n_chrs==-1)
//     {
//         fprintf(stderr, "[ERROR] -1 <hap1.p_ctg.gfa>, -2 <hap2.p_ctg.gfa>, and -n <chr_number> are required for UTG-CTG mapping and phasing.\n");
//         return 1;
//     }
//     printf("start main\n");
//     asg_t *graph = gfa_read(gfa_filename);
//     map<uint32_t, map<uint32_t, set<uint32_t>>> *bubble_chain_graph = nullptr;
    
//     uint32_t **connections_foward;
//     CALLOC(connections_foward, graph->n_seq);
//     for (int i = 0; i < graph->n_seq; i++)
//     {
//         CALLOC(connections_foward[i], graph->n_seq);
//         memset(connections_foward[i], 0, sizeof(*connections_foward[i]));
//     }
//     uint32_t **connections_backward;
//     CALLOC(connections_backward, graph->n_seq);
//     for (int i = 0; i < graph->n_seq; i++)
//     {
//         CALLOC(connections_backward[i], graph->n_seq);
//         memset(connections_backward[i], 0, sizeof(*connections_backward[i]));
//     }
//     ifstream infile(connectionFile);
//     uint32_t i, j, count_forward, count_backward;
//     while (infile >> i >> j >> count_backward >> count_forward)
//     {
//         connections_backward[i][j] = count_backward;
//         connections_backward[j][i] = count_backward;
//         connections_foward[i][j] = count_forward;
//         connections_foward[j][i] = count_forward;
//     }
    
//     std::string utg_gfa = std::string(gfa_filename);

//     // if (g_params.is_plant)
//     // {
//     //     printf("[INFO] Plant mode enabled. Using alternative phasing functions.\n");
//     //     bubble_chain_graph = phasing_plant_version(graph, string(output_directory), connections_foward, connections_backward);
//     // }
//     // else
//     // {
//         printf("[INFO] Default mode enabled. Using standard phasing functions.\n");
//         bubble_chain_graph = phasing_10_7(graph, string(output_directory), connections_foward, connections_backward);

        
//         if (g_params.debug_mode) {
//             printf("[INFO] Debug mode enabled. Executing get_haplotype_path_test_code...\n");
//             get_haplotype_path_test_code(connections_foward, connections_backward, graph, bubble_chain_graph,
//                                          output_directory, named_bubble_contigs, gfa_filename, g_params);
//         } else {
//             printf("[INFO] Executing standard model get_haplotype_path_now...\n");
//             get_haplotype_path_now(connections_foward, connections_backward, graph, bubble_chain_graph,
//                                    output_directory, named_bubble_contigs, gfa_filename, g_params);
//         }
//     // }
//     return 0;
// }




// int mapping(int argc, char *argv[])
// {
// 	return main_poreC_map_test(argc, argv);
// }


// int main(int argc, char *argv[])
// {
// 	extern double yak_realtime(void);
// 	extern double yak_cputime(void);
// 	extern void yak_reset_realtime(void);
// 	double t_start;
// 	int ret = 0, i;

// 	if (argc == 1)
// 	{
// 		fprintf(stderr, "Usage: HapFold <command> <arguments> <inputs>\n");
// 		fprintf(stderr, "Commands:\n");
// 		fprintf(stderr, "  scaffolding    		 use Hi-C/Pore-C data to resolve haplotypes\n");
// 		fprintf(stderr, "  mapping           	 map Hi-C/Pore-C data to sequences in the graph\n");
// 		fprintf(stderr, "  version               print version number\n");
// 		return 1;
// 	}
// 	yak_reset_realtime();
// 	t_start = yak_realtime();
// 	if (strcmp(argv[1], "scaffolding") == 0)
// 		ret = main_phasing_scaffolding(argc - 1, argv + 1);
// 	else if (strcmp(argv[1], "mapping") == 0)
// 		ret = mapping(argc - 1, argv + 1);
// 	else if (strcmp(argv[1], "version") == 0)
// 	{
// 		printf("HapFold: %s\n", HapFold_VERSION);
// 		return 0;
// 	}
// 	else
// 	{
// 		fprintf(stderr, "[E::%s] unknown command\n", __func__);
// 		return 1;
// 	}
// 	if (ret == 0)
// 	{
// 		fprintf(stderr, "[M::%s] Version: %s\n", __func__, HapFold_VERSION);
// 		fprintf(stderr, "[M::%s] CMD:", __func__);
// 		for (i = 0; i < argc; ++i)
// 			fprintf(stderr, " %s", argv[i]);
// 		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, yak_realtime() - t_start, yak_cputime());
// 	}
// 	return ret;
// }



/*
 * HapFold - A graph-based haplotype reconstruction framework
 */

extern "C"
{
#include "paf.h"
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "yak-priv.h"
#include "graph.h"
#include "graph_refining.h"
#include "phasing2scaffolding.h"
#include "mapping.h"

/*
 * 注意这里使用完整相对路径，避免把整个 hifiasm 目录
 * 加入 HapFold 的全局 include path。
 */
#include "hifiasm/hifiasm_entry.h"

#define HapFold_VERSION "1.0.0"

namespace hapfold
{

typedef struct
{
    int num_reads;
    char **read_file_names;
    char *output_file_name;
    int thread_num;
} ps_opt_t;

/*
 * 放在 hapfold namespace 以后，不再与 hifiasm 的 asm_opt 冲突。
 */
static ps_opt_t asm_opt;

static void init_opt(ps_opt_t *opt)
{
    memset(opt, 0, sizeof(ps_opt_t));

    opt->num_reads = 0;
    opt->read_file_names = NULL;
    opt->output_file_name = NULL;
    opt->thread_num = 1;
}

static void destory_opt(ps_opt_t *opt)
{
    if (opt->read_file_names != NULL)
    {
        free(opt->read_file_names);
        opt->read_file_names = NULL;
    }
}

/*
 * 你原来的 Print_H、mm_parse_num、
 * read_named_bubble_contigs 等继续放在这里。
 */

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

static int run_mcl_self_test()
{
    char directory_template[] = "/tmp/hapfold_mcl_selftest_XXXXXX";
    char *created = mkdtemp(directory_template);
    if (created == NULL)
    {
        perror("mkdtemp");
        return 1;
    }
    const std::string output_directory(created);

    std::vector<contig_chains> chains(8);
    for (size_t i = 0; i < chains.size(); ++i)
    {
        chains[i].index = 100 + i;
        chains[i].group_id = i < 4 ? 10 : 20;
        chains[i].group_id_new = UINT32_MAX;
        chains[i].path_length = 1000000;
        chains[i].haplo_sequences = new std::string(100, "ACGT"[i % 4]);
        chains[i].contig_info_output.push_back(
            std::make_pair(std::string("synthetic_") + std::to_string(i), true));
        chains[i].utg_path_node.push_back((uint32_t)i);
        chains[i].is_paired = true;
        chains[i].other_index = 100 + (i % 2 == 0 ? i + 1 : i - 1);
    }

    const size_t endpoint_count = chains.size() * 2;
    uint32_t **contacts = (uint32_t **)calloc(endpoint_count, sizeof(uint32_t *));
    for (size_t i = 0; i < endpoint_count; ++i)
        contacts[i] = (uint32_t *)calloc(endpoint_count, sizeof(uint32_t));

    auto set_contact = [&](size_t a, int ae, size_t b, int be, uint32_t value) {
        contacts[2 * a + ae][2 * b + be] = value;
        contacts[2 * b + be][2 * a + ae] = value;
    };
    set_contact(0, 1, 2, 0, 1000);
    set_contact(1, 1, 3, 0, 950);
    set_contact(4, 1, 6, 0, 900);
    set_contact(5, 1, 7, 0, 850);
    for (size_t a = 0; a < 4; ++a)
        for (size_t b = 4; b < 8; ++b)
            set_contact(a, 0, b, 0, 1);

    GlobalParams params;
    params.n_chrs = 2;
    params.paired_global_merge = "supported";
    params.paired_merge_min_links = 100;
    params.paired_merge_min_confidence = 1.5;
    std::ofstream fasta(output_directory + "/scaffold.fa");
    bool ok = run_pair_aware_mcl_scaffolding(
        chains, contacts, NULL, output_directory, params, fasta);
    fasta.close();

    std::map<uint32_t, int> chain_cluster;
    std::ifstream clusters(output_directory + "/chromosome_clusters.tsv");
    std::string line;
    std::getline(clusters, line);
    while (std::getline(clusters, line))
    {
        std::stringstream fields(line);
        std::string id, name, component, pair_id, cluster;
        std::getline(fields, id, '\t');
        std::getline(fields, name, '\t');
        std::getline(fields, component, '\t');
        std::getline(fields, pair_id, '\t');
        std::getline(fields, cluster, '\t');
        if (!id.empty() && !cluster.empty())
            chain_cluster[(uint32_t)strtoul(id.c_str(), NULL, 10)] = atoi(cluster.c_str());
    }

    bool constraints_ok = chain_cluster.size() == chains.size();
    std::set<int> observed_clusters;
    for (const auto &entry : chain_cluster) observed_clusters.insert(entry.second);
    constraints_ok = constraints_ok && observed_clusters.size() == 2;

    std::map<uint32_t, int> chain_scaffold;
    std::ifstream paths(output_directory + "/scaffold_mcl_paths.tsv");
    std::getline(paths, line);
    while (std::getline(paths, line))
    {
        std::stringstream fields(line);
        std::string scaffold, cluster, id, name, orientation;
        std::getline(fields, scaffold, '\t');
        std::getline(fields, cluster, '\t');
        std::getline(fields, id, '\t');
        std::getline(fields, name, '\t');
        std::getline(fields, orientation, '\t');
        if (!id.empty() && !scaffold.empty())
            chain_scaffold[(uint32_t)strtoul(id.c_str(), NULL, 10)] = atoi(scaffold.c_str());
    }
    constraints_ok = constraints_ok && chain_scaffold.size() == chains.size();
    for (size_t i = 0; i < chains.size(); i += 2)
        constraints_ok = constraints_ok &&
                         chain_cluster[chains[i].index] == chain_cluster[chains[i + 1].index] &&
                         chain_scaffold[chains[i].index] != chain_scaffold[chains[i + 1].index];

    for (size_t i = 0; i < endpoint_count; ++i) free(contacts[i]);
    free(contacts);
    for (contig_chains &chain : chains) delete chain.haplo_sequences;

    if (!ok || !constraints_ok)
    {
        fprintf(stderr, "[SELFTEST::MCL] FAILED; artifacts: %s\n", output_directory.c_str());
        return 1;
    }
    fprintf(stderr, "[SELFTEST::MCL] PASSED; artifacts: %s\n", output_directory.c_str());
    return 0;
}

/*
 * 原来的函数内容保持不变。
 */
int main_phasing_scaffolding(int argc, char *argv[])
{
    ketopt_t o = KETOPT_INIT;
    int c;

    GlobalParams g_params;
    g_params.n_chrs = -1;
    static ko_longopt_t longopts[] = {
        {"hic_scaffold_threshold_ratio", ko_required_argument, 301},
        {"debug", ko_no_argument, 302}, 
        {"chain_len_thresh", ko_required_argument, 303},     // 对应 > 12M 参与迭代的阈值
        {"scaffold_len_thresh", ko_required_argument, 304},  // 对应 > 300K 直接输出的阈值
        {"global-scaffolding", ko_required_argument, 305},
        {"paired-global-merge", ko_required_argument, 306},
        {"mcl-inflation", ko_required_argument, 307},
        {"paired-merge-min-links", ko_required_argument, 308},
        {"paired-merge-min-confidence", ko_required_argument, 309},
        {"component-seed-boost", ko_required_argument, 310},
        {"split-chain-list", ko_required_argument, 311},
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
        // else if (c == 'p') 
        //     g_params.is_plant = true;
        else if (c == 'd' || c == 302) 
            g_params.debug_mode = true;
        else if (c == 301) 
            g_params.hic_scaffold_threshold_ratio = atof(o.arg);
        else if (c == 303) // 捕获 12M 阈值
            g_params.chain_len_threshold = atoi(o.arg);
        else if (c == 304) // 捕获 300K 阈值
            g_params.scaffold_len_threshold = atoi(o.arg);
        else if (c == 305)
            g_params.global_scaffolding_mode = string(o.arg);
        else if (c == 306)
            g_params.paired_global_merge = string(o.arg);
        else if (c == 307)
            g_params.mcl_inflation = atof(o.arg);
        else if (c == 308)
            g_params.paired_merge_min_links = strtoul(o.arg, NULL, 10);
        else if (c == 309)
            g_params.paired_merge_min_confidence = atof(o.arg);
        else if (c == 310)
            g_params.component_seed_boost = atof(o.arg);
        else if (c == 311)
            g_params.split_chain_list = string(o.arg);
    }

    if (argc - o.ind < 3)
    {
        fprintf(stderr, "\nUsage: HapFold scaffolding [options] <mapping.txt> <assembly.gfa> <output_dir> -1 *.hap1.p_ctg.gfa -2 *.hap2.p_ctg.gfa -n chr_number\n\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -t INT      Number of threads [%d]\n", g_params.n_threads);
        fprintf(stderr, "  -n INT      Expected number of chromosomes (e.g., 78 for chicken) [%d]\n", g_params.n_chrs);
        fprintf(stderr, "  -e STR      Restriction enzymes separated by comma (e.g., GATC,GANTC) [%s]\n", g_params.enzymes_unsplit.c_str());
        fprintf(stderr, "  -c FILE     Path to contig_hap_nodes.txt (debug for Hi-C phasing)\n");
        fprintf(stderr, "  -u FILE     Output path/name for the UTG-CTG mapping file [default: <output_dir>/utg_ctg_mappings.csv]\n");
        fprintf(stderr, "  -1 FILE     Path to haplotype 1 GFA file (*.hap1.p_ctg.gfa)\n");
        fprintf(stderr, "  -2 FILE     Path to haplotype 2 GFA file (*.hap2.p_ctg.gfa)\n");
        fprintf(stderr, "  -i BOOL     Enable identity check on contigs (true/false) [%s]\n", (g_params.check_identity ? "true" : "false"));
        fprintf(stderr, "  -f FILE     Precomputed identity file path; if omitted, check will run automatically [%s]\n", g_params.identityFile.c_str());
        // fprintf(stderr, "  -p          Enable plant mode (uses alternative phasing algorithms) [optional]\n"); 
        fprintf(stderr, "  -d, --debug Enable debug mode to run test code functions [optional]\n"); 
        fprintf(stderr, "  --hic_scaffold_threshold_ratio FLOAT  Threshold ratio for Hi-C scaffolding [%.2f]\n", g_params.hic_scaffold_threshold_ratio);
        fprintf(stderr, "  --chain_len_thresh INT                Length threshold to join contig_chain for iterative merging [%d]\n", g_params.chain_len_threshold);
        fprintf(stderr, "  --scaffold_len_thresh INT             Length threshold to directly output to scaffold.fa [%d]\n", g_params.scaffold_len_threshold);
        fprintf(stderr, "  --global-scaffolding STR              Global method: mcl or legacy [%s]\n", g_params.global_scaffolding_mode.c_str());
        fprintf(stderr, "  --paired-global-merge STR             Pair coupling: off, supported or inferred [%s]\n", g_params.paired_global_merge.c_str());
        fprintf(stderr, "  --mcl-inflation FLOAT                 MCL inflation; 0 scans automatically [%.2f]\n", g_params.mcl_inflation);
        fprintf(stderr, "  --paired-merge-min-links INT          Minimum normalized links for a coupled edge [%u]\n", g_params.paired_merge_min_links);
        fprintf(stderr, "  --paired-merge-min-confidence FLOAT   Best/second-best endpoint ratio [%.2f]\n", g_params.paired_merge_min_confidence);
        fprintf(stderr, "  --component-seed-boost FLOAT          Soft within-component multiplier [%.2f]\n", g_params.component_seed_boost);
        fprintf(stderr, "  --split-chain-list FILE               First-contig or internal chain IDs to restore as source contigs\n");
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

    if (g_params.utg_ctg_file.empty())
    {
        std::string out_dir = std::string(output_directory);
        if (!out_dir.empty() && out_dir.back() != '/')
            out_dir += "/";

        g_params.utg_ctg_file = out_dir + "utg_ctg_mappings.csv";

        fprintf(stderr, "[INFO] No -u provided. Using default UTG-CTG mapping file: %s\n",
                g_params.utg_ctg_file.c_str());
    }
    else
    {
        fprintf(stderr, "[INFO] Using user-specified UTG-CTG mapping file: %s\n",
                g_params.utg_ctg_file.c_str());
    }


    if (g_params.hap1_gfa.empty() || g_params.hap2_gfa.empty() || g_params.n_chrs==-1)
    {
        fprintf(stderr, "[ERROR] -1 <hap1.p_ctg.gfa>, -2 <hap2.p_ctg.gfa>, and -n <chr_number> are required for UTG-CTG mapping and phasing.\n");
        return 1;
    }
    if (g_params.global_scaffolding_mode != "mcl" && g_params.global_scaffolding_mode != "legacy")
    {
        fprintf(stderr, "[ERROR] --global-scaffolding must be mcl or legacy.\n");
        return 1;
    }
    if (g_params.paired_global_merge != "off" &&
        g_params.paired_global_merge != "supported" &&
        g_params.paired_global_merge != "inferred")
    {
        fprintf(stderr, "[ERROR] --paired-global-merge must be off, supported or inferred.\n");
        return 1;
    }
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

    // if (g_params.is_plant)
    // {
    //     printf("[INFO] Plant mode enabled. Using alternative phasing functions.\n");
    //     bubble_chain_graph = phasing_plant_version(graph, string(output_directory), connections_foward, connections_backward);
    // }
    // else
    // {
        printf("[INFO] Default mode enabled. Using standard phasing functions.\n");
        bubble_chain_graph = phasing_10_7(graph, string(output_directory), connections_foward, connections_backward);

        
        if (g_params.debug_mode) {
            printf("[INFO] Debug mode enabled. Executing get_haplotype_path_test_code...\n");
            get_haplotype_path_test_code(connections_foward, connections_backward, graph, bubble_chain_graph,
                                         output_directory, named_bubble_contigs, gfa_filename, g_params);
        } else {
            printf("[INFO] Executing standard model get_haplotype_path_now...\n");
            get_haplotype_path_now(connections_foward, connections_backward, graph, bubble_chain_graph,
                                   output_directory, named_bubble_contigs, gfa_filename, g_params);
        }
    // }
    return 0;
}

int mapping_entry(int argc, char *argv[])
{
    return main_poreC_map_test(argc, argv);
}

static bool path_exists(const std::string &path)
{
    struct stat st;
    return !path.empty() && stat(path.c_str(), &st) == 0;
}

static std::string parent_path(const std::string &path)
{
    const std::string::size_type pos = path.find_last_of('/');
    if (pos == std::string::npos)
        return "";
    if (pos == 0)
        return "/";
    return path.substr(0, pos);
}

static bool make_directories(const std::string &path)
{
    if (path.empty() || path == "." || path_exists(path))
        return true;

    const std::string parent = parent_path(path);
    if (!parent.empty() && parent != path && !make_directories(parent))
        return false;

    if (mkdir(path.c_str(), 0775) == 0 || errno == EEXIST)
        return true;

    fprintf(stderr, "[ERROR] Cannot create directory %s: %s\n",
            path.c_str(), strerror(errno));
    return false;
}

static bool create_owned_work_directory(const std::string &path)
{
    const std::string marker = path + "/.hapfold-owned-workdir";
    if (path_exists(path) && !path_exists(marker))
    {
        fprintf(stderr,
                "[ERROR] Refusing to use pre-existing unowned work directory: %s\n",
                path.c_str());
        return false;
    }
    if (!make_directories(path))
        return false;

    std::ofstream marker_output(marker.c_str(), std::ios::out | std::ios::trunc);
    if (!marker_output)
    {
        fprintf(stderr, "[ERROR] Cannot create work-directory marker: %s\n",
                marker.c_str());
        return false;
    }
    marker_output << "HapFold managed temporary directory\n";
    return static_cast<bool>(marker_output);
}

static bool remove_owned_tree_impl(const std::string &path)
{
    struct stat st;
    if (lstat(path.c_str(), &st) != 0)
        return errno == ENOENT;

    if (!S_ISDIR(st.st_mode) || S_ISLNK(st.st_mode))
    {
        if (unlink(path.c_str()) == 0 || errno == ENOENT)
            return true;
        fprintf(stderr, "[ERROR] Cannot remove temporary file %s: %s\n",
                path.c_str(), strerror(errno));
        return false;
    }

    DIR *directory = opendir(path.c_str());
    if (directory == NULL)
    {
        fprintf(stderr, "[ERROR] Cannot open temporary directory %s: %s\n",
                path.c_str(), strerror(errno));
        return false;
    }

    bool ok = true;
    struct dirent *entry;
    while ((entry = readdir(directory)) != NULL)
    {
        if (strcmp(entry->d_name, ".") == 0 ||
            strcmp(entry->d_name, "..") == 0)
            continue;
        if (!remove_owned_tree_impl(path + "/" + entry->d_name))
            ok = false;
    }
    closedir(directory);

    if (ok && rmdir(path.c_str()) != 0 && errno != ENOENT)
    {
        fprintf(stderr, "[ERROR] Cannot remove temporary directory %s: %s\n",
                path.c_str(), strerror(errno));
        ok = false;
    }
    return ok;
}

static bool remove_owned_work_directory(const std::string &path)
{
    const std::string marker = path + "/.hapfold-owned-workdir";
    if (path.empty() || path == "/" || path == "." || !path_exists(marker))
    {
        fprintf(stderr,
                "[ERROR] Refusing to clean unmarked work directory: %s\n",
                path.c_str());
        return false;
    }
    return remove_owned_tree_impl(path);
}

static bool require_file(const std::string &path, const char *description)
{
    if (path_exists(path))
        return true;
    fprintf(stderr, "[ERROR] Missing %s: %s\n", description, path.c_str());
    return false;
}

static bool gfa_to_fasta(const std::string &gfa_path,
                         const std::string &fasta_path)
{
    std::ifstream input(gfa_path.c_str());
    if (!input)
    {
        fprintf(stderr, "[ERROR] Cannot open unitig GFA: %s\n", gfa_path.c_str());
        return false;
    }

    std::ofstream output(fasta_path.c_str());
    if (!output)
    {
        fprintf(stderr, "[ERROR] Cannot create unitig FASTA: %s\n",
                fasta_path.c_str());
        return false;
    }

    std::string line;
    size_t sequence_count = 0;
    size_t missing_sequence_count = 0;
    while (std::getline(input, line))
    {
        if (line.size() < 2 || line[0] != 'S' || line[1] != '\t')
            continue;

        std::string name, sequence;
        std::stringstream fields(line);
        std::string record_type;
        std::getline(fields, record_type, '\t');
        std::getline(fields, name, '\t');
        std::getline(fields, sequence, '\t');

        if (name.empty() || sequence.empty() || sequence == "*")
        {
            ++missing_sequence_count;
            continue;
        }

        output << '>' << name << '\n' << sequence << '\n';
        ++sequence_count;
    }

    if (!output || sequence_count == 0)
    {
        fprintf(stderr,
                "[ERROR] No segment sequences were extracted from %s\n",
                gfa_path.c_str());
        return false;
    }

    fprintf(stderr,
            "[M::run_pipeline] Extracted %zu unitig sequences to %s",
            sequence_count, fasta_path.c_str());
    if (missing_sequence_count != 0)
        fprintf(stderr, " (%zu sequence-less segments skipped)",
                missing_sequence_count);
    fputc('\n', stderr);
    return true;
}

static void print_run_help()
{
    fprintf(stderr,
            "\nUsage:\n"
            "  HapFold run -1 HIC_R1 -2 HIC_R2 -n CHROMOSOMES [options]\n"
            "              -- [native hifiasm options] <assembly_reads>\n\n"
            "Required run options:\n"
            "  -1, --hic1 FILE          Hi-C/Pore-C read 1\n"
            "  -2, --hic2 FILE          Hi-C/Pore-C read 2\n"
            "  -n, --chromosomes INT    Expected chromosome count\n\n"
            "Pipeline options:\n"
            "  -t, --threads INT        Threads shared by all three stages [32]\n"
            "  -o, --output-prefix STR  hifiasm/HapFold shared prefix [hifiasm.asm]\n"
            "      --hapfold-output DIR  Scaffolding directory [<prefix>.hapfold]\n"
            "      --hifiasm-mode STR     Hifiasm mode: default, trio, or hic [hic]\n"
            "      --hifiasm-hap1-yak FILE  Paternal/haplotype-1 yak dump (trio mode)\n"
            "      --hifiasm-hap2-yak FILE  Maternal/haplotype-2 yak dump (trio mode)\n"
            "      --keep-hifiasm-output Keep raw Hifiasm files and restart caches\n"
            "      --high-quality-utg    Recompute with 5 correction rounds and enhanced cleaning\n"
            "      --mapping-output FILE  Mapping result [DIR/mapping.txt]\n"
            "      --utg-gfa FILE         Override auto-detected p_utg.gfa\n"
            "      --hap1-gfa FILE        Override auto-detected hap1.p_ctg.gfa\n"
            "      --hap2-gfa FILE        Override auto-detected hap2.p_ctg.gfa\n"
            "      --utg-fasta FILE       Unitig FASTA to create/use [DIR/p_utg.fa]\n"
            "  -e STR                    Restriction enzymes for scaffolding\n"
            "  -i BOOL                   Scaffolding identity check\n"
            "  -f FILE                   Precomputed identity file\n"
            "  -c FILE                   contig_hap_nodes file\n"
            "  -u FILE                   UTG-CTG mapping output path\n"
            "  -d, --debug               Scaffolding debug mode\n"
            "      --hic_scaffold_threshold_ratio FLOAT\n"
            "      --chain_len_thresh INT\n"
            "      --scaffold_len_thresh INT\n"
            "      --global-scaffolding mcl|legacy [mcl]\n"
            "      --paired-global-merge off|supported|inferred [supported]\n"
            "      --mcl-inflation FLOAT [automatic scan]\n"
            "      --paired-merge-min-links INT [100]\n"
            "      --paired-merge-min-confidence FLOAT [1.5]\n"
            "      --component-seed-boost FLOAT [1.2]\n\n"
            "Example:\n"
            "  HapFold run -1 hic.R1.fq.gz -2 hic.R2.fq.gz -n 46 -t 32 \\\n"
            "    -o result/asm --high-quality-utg -- hifi.fq.gz\n\n");
}

int run_pipeline(int argc, char *argv[])
{
    int separator_index = -1;

    /*
     * 此时：
     *
     * argv[0] = "run"
     */
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "--") == 0)
        {
            separator_index = i;
            break;
        }
    }

    if (separator_index < 0)
    {
        fprintf(stderr, "[ERROR] run requires '--' before hifiasm arguments.\n");
        print_run_help();
        return 1;
    }

    if (separator_index + 1 >= argc)
    {
        fprintf(stderr,
                "[ERROR] No hifiasm arguments were supplied "
                "after '--'.\n");

        return 1;
    }

    std::string hic1, hic2, output_prefix = "hifiasm.asm";
    std::string output_dir, mapping_output;
    std::string utg_gfa, hap1_gfa, hap2_gfa, utg_fasta;
    std::string enzymes, identity_check, identity_file;
    std::string contig_hap_file, utg_ctg_file;
    std::string hic_threshold, chain_threshold, scaffold_threshold;
    std::string global_scaffolding, paired_global_merge, mcl_inflation;
    std::string paired_merge_min_links, paired_merge_min_confidence, component_seed_boost;
    std::string hifiasm_mode = "hic";
    std::string hifiasm_hap1_yak, hifiasm_hap2_yak;
    int threads = 32;
    int chromosomes = -1;
    bool debug_mode = false;
    bool high_quality_utg = false;
    bool keep_hifiasm_output = false;

    for (int i = 1; i < separator_index; ++i)
    {
        const std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help")
        {
            print_run_help();
            return 0;
        }
        if (arg == "-d" || arg == "--debug")
        {
            debug_mode = true;
            continue;
        }
        if (arg == "--high-quality-utg")
        {
            high_quality_utg = true;
            continue;
        }
        if (arg == "--keep-hifiasm-output")
        {
            keep_hifiasm_output = true;
            continue;
        }

        if (i + 1 >= separator_index)
        {
            fprintf(stderr, "[ERROR] Missing value for run option %s\n", argv[i]);
            return 1;
        }
        const std::string value(argv[++i]);

        if (arg == "-1" || arg == "--hic1") hic1 = value;
        else if (arg == "-2" || arg == "--hic2") hic2 = value;
        else if (arg == "-n" || arg == "--chromosomes") chromosomes = atoi(value.c_str());
        else if (arg == "-t" || arg == "--threads") threads = atoi(value.c_str());
        else if (arg == "-o" || arg == "--output-prefix") output_prefix = value;
        else if (arg == "--hapfold-output") output_dir = value;
        else if (arg == "--hifiasm-mode") hifiasm_mode = value;
        else if (arg == "--hifiasm-hap1-yak") hifiasm_hap1_yak = value;
        else if (arg == "--hifiasm-hap2-yak") hifiasm_hap2_yak = value;
        else if (arg == "--mapping-output") mapping_output = value;
        else if (arg == "--utg-gfa") utg_gfa = value;
        else if (arg == "--hap1-gfa") hap1_gfa = value;
        else if (arg == "--hap2-gfa") hap2_gfa = value;
        else if (arg == "--utg-fasta") utg_fasta = value;
        else if (arg == "-e") enzymes = value;
        else if (arg == "-i") identity_check = value;
        else if (arg == "-f") identity_file = value;
        else if (arg == "-c") contig_hap_file = value;
        else if (arg == "-u") utg_ctg_file = value;
        else if (arg == "--hic_scaffold_threshold_ratio") hic_threshold = value;
        else if (arg == "--chain_len_thresh") chain_threshold = value;
        else if (arg == "--scaffold_len_thresh") scaffold_threshold = value;
        else if (arg == "--global-scaffolding") global_scaffolding = value;
        else if (arg == "--paired-global-merge") paired_global_merge = value;
        else if (arg == "--mcl-inflation") mcl_inflation = value;
        else if (arg == "--paired-merge-min-links") paired_merge_min_links = value;
        else if (arg == "--paired-merge-min-confidence") paired_merge_min_confidence = value;
        else if (arg == "--component-seed-boost") component_seed_boost = value;
        else
        {
            fprintf(stderr, "[ERROR] Unknown run option: %s\n", arg.c_str());
            print_run_help();
            return 1;
        }
    }

    if (hic1.empty() || hic2.empty() || chromosomes <= 0 || threads <= 0)
    {
        fprintf(stderr,
                "[ERROR] run requires -1, -2, and a positive -n; "
                "-t must also be positive.\n");
        print_run_help();
        return 1;
    }
    if (!require_file(hic1, "Hi-C read 1") ||
        !require_file(hic2, "Hi-C read 2"))
        return 1;
    if (hifiasm_mode != "default" &&
        hifiasm_mode != "trio" &&
        hifiasm_mode != "hic")
    {
        fprintf(stderr,
                "[ERROR] --hifiasm-mode must be default, trio, or hic.\n");
        return 1;
    }
    if (hifiasm_mode == "trio")
    {
        if (hifiasm_hap1_yak.empty() || hifiasm_hap2_yak.empty())
        {
            fprintf(stderr,
                    "[ERROR] Trio mode requires --hifiasm-hap1-yak and "
                    "--hifiasm-hap2-yak.\n");
            return 1;
        }
        if (!require_file(hifiasm_hap1_yak, "haplotype-1 yak dump") ||
            !require_file(hifiasm_hap2_yak, "haplotype-2 yak dump"))
            return 1;
    }
    else if (!hifiasm_hap1_yak.empty() || !hifiasm_hap2_yak.empty())
    {
        fprintf(stderr,
                "[ERROR] --hifiasm-hap1-yak/--hifiasm-hap2-yak are only "
                "valid with --hifiasm-mode trio.\n");
        return 1;
    }

    std::vector<std::string> native_hifiasm_args;
    for (int i = separator_index + 1; i < argc; ++i)
    {
        const std::string arg(argv[i]);
        if (arg == "-o" || arg == "-t" ||
            (arg.size() > 2 &&
             (arg.substr(0, 2) == "-o" || arg.substr(0, 2) == "-t")))
        {
            fprintf(stderr,
                    "[ERROR] %s is a shared run option. Put -o/-t before "
                    "'--' and specify each only once.\n",
                    arg.c_str());
            return 1;
        }
        if (arg == "--h1" || arg == "--h2" ||
            arg == "-1" || arg == "-2" || arg == "-3" || arg == "-4")
        {
            fprintf(stderr,
                    "[ERROR] Hifiasm phasing inputs are managed by "
                    "--hifiasm-mode; do not put %s after '--'.\n",
                    arg.c_str());
            return 1;
        }
        native_hifiasm_args.push_back(arg);
    }

    if (output_dir.empty())
        output_dir = output_prefix + ".hapfold";
    /*
     * Scaffolding rebuilds output_dir from scratch (phasing_10_7 starts by
     * removing that directory).  Keeping the embedded-hifiasm products below
     * output_dir therefore deleted the GFA and mapping inputs immediately
     * before scaffolding tried to open them.  Use a managed sibling directory
     * so scaffolding can safely recreate output_dir without invalidating its
     * own inputs.
     */
    const std::string work_dir = output_dir + ".hapfold-work";
    const std::string hifiasm_work_dir =
        work_dir + "/hifiasm-" + hifiasm_mode;
    const std::string hifiasm_output_prefix =
        keep_hifiasm_output ? output_prefix : hifiasm_work_dir + "/asm";

    if (mapping_output.empty())
        mapping_output = output_dir + "/mapping.txt";
    if (utg_fasta.empty())
        utg_fasta = work_dir + "/p_utg.fa";

    if (!make_directories(parent_path(output_prefix)) ||
        !make_directories(output_dir) ||
        !create_owned_work_directory(work_dir) ||
        (!keep_hifiasm_output && !make_directories(hifiasm_work_dir)) ||
        !make_directories(parent_path(mapping_output)) ||
        !make_directories(parent_path(utg_fasta)))
        return 1;

    const auto report_preserved_work = [&]() {
        fprintf(stderr,
                "[M::run_pipeline] Preserving managed work directory after "
                "failure: %s\n",
                work_dir.c_str());
    };

    std::vector<std::string> hifiasm_storage;
    hifiasm_storage.push_back("hifiasm");
    hifiasm_storage.push_back("-o");
    hifiasm_storage.push_back(hifiasm_output_prefix);
    hifiasm_storage.push_back("-t");
    hifiasm_storage.push_back(std::to_string(threads));
    if (high_quality_utg)
    {
        /*
         * Conservative quality preset: increase correction, graph-cleaning
         * iterations, and the fallback overlap candidates. Native hifiasm
         * arguments are appended later and may override -r/-a/-N.
         */
        hifiasm_storage.push_back("-r");
        hifiasm_storage.push_back("5");
        hifiasm_storage.push_back("-a");
        hifiasm_storage.push_back("6");
        hifiasm_storage.push_back("-N");
        hifiasm_storage.push_back("150");
        hifiasm_storage.push_back("-i");
    }
    hifiasm_storage.insert(hifiasm_storage.end(),
                            native_hifiasm_args.begin(),
                            native_hifiasm_args.end());
    if (hifiasm_mode == "hic")
    {
        hifiasm_storage.push_back("--h1");
        hifiasm_storage.push_back(hic1);
        hifiasm_storage.push_back("--h2");
        hifiasm_storage.push_back(hic2);
    }
    else if (hifiasm_mode == "trio")
    {
        hifiasm_storage.push_back("-1");
        hifiasm_storage.push_back(hifiasm_hap1_yak);
        hifiasm_storage.push_back("-2");
        hifiasm_storage.push_back(hifiasm_hap2_yak);
    }
    if (!keep_hifiasm_output)
        hifiasm_storage.push_back("--no-write-cache");

    std::vector<char *> hifiasm_argv;
    for (size_t i = 0; i < hifiasm_storage.size(); ++i)
        hifiasm_argv.push_back(const_cast<char *>(hifiasm_storage[i].c_str()));
    hifiasm_argv.push_back(nullptr);

    const int hifiasm_argc =
        static_cast<int>(hifiasm_argv.size()) - 1;

    fprintf(stderr,
            "[M::run_pipeline] Starting embedded hifiasm "
            "(mode: %s, prefix: %s, raw output: %s)\n",
            hifiasm_mode.c_str(), hifiasm_output_prefix.c_str(),
            keep_hifiasm_output ? "kept" : "temporary");
    if (high_quality_utg)
        fprintf(stderr,
                "[M::run_pipeline] High-quality UTG preset: "
                "-r 5 -a 6 -N 150 -i\n");

    int ret = ::hifiasm_main(
        hifiasm_argc,
        hifiasm_argv.data());

    if (ret != 0)
    {
        fprintf(stderr,
                "[ERROR] Embedded hifiasm failed "
                "with exit code %d\n",
                ret);
        report_preserved_work();
        return ret;
    }

    fprintf(stderr,
            "[M::run_pipeline] Embedded hifiasm completed\n");

    const std::string mode_suffix =
        hifiasm_mode == "default" ? "bp" :
        hifiasm_mode == "trio" ? "dip" : "hic";
    if (utg_gfa.empty())
        utg_gfa = hifiasm_output_prefix + "." + mode_suffix + ".p_utg.gfa";
    if (hap1_gfa.empty())
        hap1_gfa = hifiasm_output_prefix + "." + mode_suffix +
                   ".hap1.p_ctg.gfa";
    if (hap2_gfa.empty())
        hap2_gfa = hifiasm_output_prefix + "." + mode_suffix +
                   ".hap2.p_ctg.gfa";

    if (!require_file(utg_gfa, "hifiasm unitig GFA") ||
        !require_file(hap1_gfa, "hifiasm haplotype-1 GFA") ||
        !require_file(hap2_gfa, "hifiasm haplotype-2 GFA"))
    {
        fprintf(stderr,
                "[ERROR] Use --utg-gfa/--hap1-gfa/--hap2-gfa when "
                "hifiasm produced non-standard names.\n");
        report_preserved_work();
        return 1;
    }

    if (!gfa_to_fasta(utg_gfa, utg_fasta))
    {
        report_preserved_work();
        return 1;
    }

    const std::string threads_string = std::to_string(threads);
    std::vector<std::string> mapping_storage;
    mapping_storage.push_back("mapping");
    mapping_storage.push_back("-t");
    mapping_storage.push_back(threads_string);
    mapping_storage.push_back("-1");
    mapping_storage.push_back(hic1);
    mapping_storage.push_back("-2");
    mapping_storage.push_back(hic2);
    mapping_storage.push_back("-o");
    mapping_storage.push_back(mapping_output);
    mapping_storage.push_back(utg_fasta);

    std::vector<char *> mapping_argv;
    for (size_t i = 0; i < mapping_storage.size(); ++i)
        mapping_argv.push_back(const_cast<char *>(mapping_storage[i].c_str()));
    mapping_argv.push_back(nullptr);

    fprintf(stderr, "[M::run_pipeline] Starting mapping -> %s\n",
            mapping_output.c_str());
    ret = mapping_entry(static_cast<int>(mapping_storage.size()),
                        mapping_argv.data());
    if (ret != 0)
    {
        fprintf(stderr, "[ERROR] Mapping failed with exit code %d\n", ret);
        report_preserved_work();
        return ret;
    }
    if (!require_file(mapping_output, "mapping output"))
    {
        report_preserved_work();
        return 1;
    }
    fprintf(stderr, "[M::run_pipeline] Mapping completed\n");

    std::vector<std::string> scaffold_storage;
    scaffold_storage.push_back("scaffolding");
    scaffold_storage.push_back("-t");
    scaffold_storage.push_back(threads_string);
    scaffold_storage.push_back("-n");
    scaffold_storage.push_back(std::to_string(chromosomes));
    scaffold_storage.push_back("-1");
    scaffold_storage.push_back(hap1_gfa);
    scaffold_storage.push_back("-2");
    scaffold_storage.push_back(hap2_gfa);
    if (!enzymes.empty()) {
        scaffold_storage.push_back("-e"); scaffold_storage.push_back(enzymes);
    }
    if (!identity_check.empty()) {
        scaffold_storage.push_back("-i"); scaffold_storage.push_back(identity_check);
    }
    if (!identity_file.empty()) {
        scaffold_storage.push_back("-f"); scaffold_storage.push_back(identity_file);
    }
    if (!contig_hap_file.empty()) {
        scaffold_storage.push_back("-c"); scaffold_storage.push_back(contig_hap_file);
    }
    if (!utg_ctg_file.empty()) {
        scaffold_storage.push_back("-u"); scaffold_storage.push_back(utg_ctg_file);
    }
    if (debug_mode)
        scaffold_storage.push_back("-d");
    if (!hic_threshold.empty()) {
        scaffold_storage.push_back("--hic_scaffold_threshold_ratio");
        scaffold_storage.push_back(hic_threshold);
    }
    if (!chain_threshold.empty()) {
        scaffold_storage.push_back("--chain_len_thresh");
        scaffold_storage.push_back(chain_threshold);
    }
    if (!scaffold_threshold.empty()) {
        scaffold_storage.push_back("--scaffold_len_thresh");
        scaffold_storage.push_back(scaffold_threshold);
    }
    if (!global_scaffolding.empty()) {
        scaffold_storage.push_back("--global-scaffolding");
        scaffold_storage.push_back(global_scaffolding);
    }
    if (!paired_global_merge.empty()) {
        scaffold_storage.push_back("--paired-global-merge");
        scaffold_storage.push_back(paired_global_merge);
    }
    if (!mcl_inflation.empty()) {
        scaffold_storage.push_back("--mcl-inflation");
        scaffold_storage.push_back(mcl_inflation);
    }
    if (!paired_merge_min_links.empty()) {
        scaffold_storage.push_back("--paired-merge-min-links");
        scaffold_storage.push_back(paired_merge_min_links);
    }
    if (!paired_merge_min_confidence.empty()) {
        scaffold_storage.push_back("--paired-merge-min-confidence");
        scaffold_storage.push_back(paired_merge_min_confidence);
    }
    if (!component_seed_boost.empty()) {
        scaffold_storage.push_back("--component-seed-boost");
        scaffold_storage.push_back(component_seed_boost);
    }
    scaffold_storage.push_back(mapping_output);
    scaffold_storage.push_back(utg_gfa);
    scaffold_storage.push_back(output_dir);

    std::vector<char *> scaffold_argv;
    for (size_t i = 0; i < scaffold_storage.size(); ++i)
        scaffold_argv.push_back(const_cast<char *>(scaffold_storage[i].c_str()));
    scaffold_argv.push_back(nullptr);

    fprintf(stderr, "[M::run_pipeline] Starting scaffolding -> %s\n",
            output_dir.c_str());
    ret = main_phasing_scaffolding(
        static_cast<int>(scaffold_storage.size()), scaffold_argv.data());
    if (ret != 0)
    {
        fprintf(stderr, "[ERROR] Scaffolding failed with exit code %d\n", ret);
        report_preserved_work();
        return ret;
    }
    fprintf(stderr, "[M::run_pipeline] Pipeline completed successfully\n");

    if (!remove_owned_work_directory(work_dir))
    {
        fprintf(stderr,
                "[ERROR] Pipeline outputs are complete, but temporary "
                "work-directory cleanup failed: %s\n",
                work_dir.c_str());
        return 1;
    }
    fprintf(stderr, "[M::run_pipeline] Removed temporary work directory: %s\n",
            work_dir.c_str());

    return 0;
}

static void print_main_help()
{
    fprintf(stderr,
            "\n"
            "Usage: HapFold <command> [arguments]\n\n"
            "Commands:\n"
            "  hifiasm       Run the embedded hifiasm assembler\n"
            "  mapping       Map Hi-C/Pore-C reads to graph sequences\n"
            "  scaffolding   Refine graph, phase and scaffold\n"
            "  run           Run hifiasm followed by the HapFold workflow\n"
            "  selftest-mcl  Run the synthetic pair-aware MCL regression test\n"
            "  version       Print version information\n\n"
            "Examples:\n"
            "  HapFold hifiasm -o asm -t 64 reads.fastq.gz\n"
            "  HapFold mapping ...\n"
            "  HapFold scaffolding ...\n"
            "  HapFold run -1 hic.R1.fq.gz -2 hic.R2.fq.gz -n 46 "
            "-t 64 -o asm --high-quality-utg -- reads.fastq.gz\n\n");
}

} // namespace hapfold


/*
 * 整个项目中唯一真正的 main。
 */
int main(int argc, char *argv[])
{
    extern double yak_realtime(void);
    extern double yak_cputime(void);
    extern void yak_reset_realtime(void);

    int ret = 0;
    int i = 0;

    if (argc == 1)
    {
        hapfold::print_main_help();
        return 1;
    }

    /*
     * hifiasm 自己会完成 timer reset 和日志输出。
     * 因此这里直接返回，不再进入 HapFold 的通用计时部分。
     */
    if (strcmp(argv[1], "hifiasm") == 0)
    {
        return ::hifiasm_main(argc - 1, argv + 1);
    }

    /*
     * 第一版 run 内部先调用一次 hifiasm。
     */
    if (strcmp(argv[1], "run") == 0)
    {
        return hapfold::run_pipeline(argc - 1, argv + 1);
    }
    if (strcmp(argv[1], "selftest-mcl") == 0)
    {
        return hapfold::run_mcl_self_test();
    }

    yak_reset_realtime();
    const double t_start = yak_realtime();

    if (strcmp(argv[1], "scaffolding") == 0)
    {
        ret = hapfold::main_phasing_scaffolding(
            argc - 1,
            argv + 1);
    }
    else if (strcmp(argv[1], "mapping") == 0)
    {
        ret = hapfold::mapping_entry(
            argc - 1,
            argv + 1);
    }
    else if (strcmp(argv[1], "version") == 0)
    {
        printf("HapFold: %s\n", HapFold_VERSION);
        return 0;
    }
    else
    {
        fprintf(stderr,
                "[E::%s] Unknown command: %s\n",
                __func__,
                argv[1]);

        hapfold::print_main_help();
        return 1;
    }

    if (ret == 0)
    {
        fprintf(stderr,
                "[M::%s] Version: %s\n",
                __func__,
                HapFold_VERSION);

        fprintf(stderr,
                "[M::%s] CMD:",
                __func__);

        for (i = 0; i < argc; ++i)
        {
            fprintf(stderr, " %s", argv[i]);
        }

        fprintf(stderr,
                "\n[M::%s] Real time: %.3f sec; "
                "CPU: %.3f sec\n",
                __func__,
                yak_realtime() - t_start,
                yak_cputime());
    }

    return ret;
}
