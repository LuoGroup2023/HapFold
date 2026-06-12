/*
 * HapFold - A graph-based haplotype reconstruction framework
 * Copyright (C) 2024 Yichen Li 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * ...
 */


extern "C"
{
#include "paf.h"
}
#include <string.h>
#include <iostream>
#include "graph.h"
#include "sys.cpp"
#include <fstream>
#include "graph_refining.h"
#include "phasing2scaffolding.h"
#include "mapping.h"
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
    fprintf(stderr, "\nUsage: HapFold <command> [options]\n\n");
    
    fprintf(stderr, "Commands:\n");
    fprintf(stderr, "    mapping       Map Hi-C/Pore-C reads to the unitig sequences\n");
    fprintf(stderr, "    scaffolding   Refine graph, phase haplotypes, and build scaffolds\n\n");
    
    fprintf(stderr, "Global Options:\n");
    fprintf(stderr, "    -o FILE       prefix of output files/directory [%s]\n", asm_opt->output_file_name);
    fprintf(stderr, "    -t INT        number of threads [%d]\n", asm_opt->thread_num);
    fprintf(stderr, "    --version     show version number\n");
    fprintf(stderr, "    -h            show help information\n\n");

    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  Step 1. Hi-C/Pore-C Mapping:\n");
    fprintf(stderr, "    ./HapFold mapping -t 32 -1 hic_1.fq.gz -2 hic_2.fq.gz -o mapping.txt utg.fa\n\n");
    
    fprintf(stderr, "  Step 2. Scaffolding (Graph refining & Phasing):\n");
    fprintf(stderr, "    ./HapFold scaffolding mapping.txt assembly.gfa out_dir -1 hap1.p_ctg.gfa -2 hap2.p_ctg.gfa -u utg_ctg.csv\n\n");
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
int main_phasing_scaffolding(int argc, char *argv[])
{
    ketopt_t o = KETOPT_INIT;
    int c;

    GlobalParams g_params;

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
    }

    if (argc - o.ind < 3)
    {
        fprintf(stderr, "\nUsage: HapFold scaffolding [options] <mapping.txt> <assembly.gfa> <output_dir> -1 *.hap1.p_ctg.gfa -2 *.hap2.p_ctg.gfa\n\n");
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


    if (g_params.hap1_gfa.empty() || g_params.hap2_gfa.empty())
    {
        fprintf(stderr, "[ERROR] Both -1 <hap1.p_ctg.gfa> and -2 <hap2.p_ctg.gfa> are required for UTG-CTG mapping and phasing.\n");
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




int mapping(int argc, char *argv[])
{
	return main_poreC_map_test(argc, argv);
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
		fprintf(stderr, "  scaffolding    		 use Hi-C/Pore-C data to resolve haplotypes\n");
		fprintf(stderr, "  mapping           	 map Hi-C/Pore-C data to sequences in the graph\n");
		fprintf(stderr, "  version               print version number\n");
		return 1;
	}
	yak_reset_realtime();
	t_start = yak_realtime();
	if (strcmp(argv[1], "scaffolding") == 0)
		ret = main_phasing_scaffolding(argc - 1, argv + 1);
	else if (strcmp(argv[1], "mapping") == 0)
		ret = mapping(argc - 1, argv + 1);
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
