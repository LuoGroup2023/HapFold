#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
calc_HPI_via_Map.py

计算单倍型纯度指数 (HPI)。
必须提供一个 Reference Mapping 文件来明确指定参考序列的染色体和单倍型归属。

Usage:
    python calc_HPI_via_Map.py -i aln.paf -m ref_map.txt -o report.tsv
"""

import sys
import argparse
import os
from collections import defaultdict

# 尝试导入 pysam
try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate Haplotype Purity Index (HPI) using an explicit reference map.")
    
    # 核心输入
    parser.add_argument("-i", "--input", required=True, help="Input alignment file (PAF, BAM, or SAM).")
    parser.add_argument("-m", "--ref_map", required=True, help="[Required] Reference mapping file. Format: RefName <tab> Chrom <tab> Hap")
    parser.add_argument("-o", "--output", default="scaffold_HPI_report.tsv", help="Output report file path.")
    
    # 权重参数
    parser.add_argument("--w_switch", type=float, default=0.5, help="Penalty weight for haplotype switch (intra-chromosomal). Default: 0.5")
    parser.add_argument("--w_chimera", type=float, default=1.0, help="Penalty weight for chimera (inter-chromosomal). Default: 1.0")
    
    # 过滤参数
    parser.add_argument("--min_aln_len", type=int, default=10000, help="Minimum alignment length (bp). Default: 10000")
    parser.add_argument("--min_mapq", type=int, default=30, help="Minimum MAPQ. Default: 30")
    
    return parser.parse_args()

def load_ref_mapping(map_file):
    """
    加载用户提供的映射文件。
    返回字典: { ref_name : (chrom, hap) }
    """
    mapping = {}
    print(f"[Info] Loading reference map from: {map_file}", file=sys.stderr)
    try:
        with open(map_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"): continue
                
                parts = line.split()
                if len(parts) < 3:
                    print(f"[Warning] Line {line_num} in map file is invalid (needs 3 cols): {line}", file=sys.stderr)
                    continue
                
                ref_name, chrom, hap = parts[0], parts[1], parts[2]
                mapping[ref_name] = (chrom, hap)
                
    except Exception as e:
        print(f"[Error] Failed to read map file: {e}", file=sys.stderr)
        sys.exit(1)
        
    print(f"[Info] Loaded {len(mapping)} reference sequences.", file=sys.stderr)
    return mapping

def read_alignments(filepath, min_len, min_mapq):
    """
    统一读取 PAF/BAM/SAM，yield (scaffold_name, ref_name, aln_len)
    """
    if filepath.endswith(".paf") or filepath.endswith(".paf.gz"):
        opener = open
        if filepath.endswith(".gz"):
            import gzip
            opener = gzip.open
        try:
            with opener(filepath, 'rt') as f:
                for line in f:
                    cols = line.strip().split('\t')
                    if len(cols) < 12: continue
                    
                    scaff = cols[0]
                    ref = cols[5]
                    mapq = int(cols[11])
                    length = int(cols[10]) # Block length
                    
                    if mapq >= min_mapq and length >= min_len:
                        yield scaff, ref, length
        except Exception as e:
             print(f"[Error] Reading PAF: {e}", file=sys.stderr)
             sys.exit(1)

    elif filepath.endswith((".bam", ".sam", ".cram")):
        if not HAS_PYSAM:
            print("[Error] Input is BAM/SAM but 'pysam' is not installed.", file=sys.stderr)
            sys.exit(1)
        try:
            save_mode = "rb" if filepath.endswith(".bam") else "r"
            with pysam.AlignmentFile(filepath, save_mode) as samfile:
                for read in samfile.fetch(until_eof=True):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
                    if read.mapping_quality < min_mapq: continue
                    if read.reference_length < min_len: continue
                    
                    yield read.query_name, read.reference_name, read.reference_length
        except Exception as e:
            print(f"[Error] Reading BAM/SAM: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("[Error] Unsupported file format.", file=sys.stderr)
        sys.exit(1)

def main():
    args = parse_args()
    
    # 1. 加载 Map
    ref_map = load_ref_mapping(args.ref_map)
    
    # 2. 读取比对并统计
    # Data: scaffold -> { (chrom, hap) : length }
    scaffold_stats = defaultdict(lambda: defaultdict(int))
    
    print(f"[Info] Parsing alignments from: {args.input}", file=sys.stderr)
    count = 0
    mapped_records = 0
    unmapped_ref_warning = set()
    
    for scaff, ref, length in read_alignments(args.input, args.min_aln_len, args.min_mapq):
        count += 1
        
        # 核心逻辑：直接查表
        if ref in ref_map:
            chrom, hap = ref_map[ref]
            scaffold_stats[scaff][(chrom, hap)] += length
            mapped_records += 1
        else:
            # 如果比对到的 ref 不在 map 文件里，记录一下（只需记录前几个防止刷屏）
            if len(unmapped_ref_warning) < 5:
                unmapped_ref_warning.add(ref)
    
    if len(unmapped_ref_warning) > 0:
        print(f"[Warning] Some reference names in alignment were NOT found in map file (e.g., {list(unmapped_ref_warning)[0]}). These alignments are ignored.", file=sys.stderr)
    
    print(f"[Info] Processed {count} valid alignment records. Used {mapped_records} records for HPI calculation.", file=sys.stderr)
    
    # 3. 计算 HPI 并输出
    try:
        out_f = open(args.output, 'w')
    except IOError as e:
        print(f"[Error] Cannot open output: {e}", file=sys.stderr)
        sys.exit(1)
    
    # 写入表头
    header = [
        "Scaffold_Name", "Dominant_ID", "Total_Aligned_Bp",
        "Pure_Bp", "Switch_Bp", "Chimera_Bp",
        "Pure_Ratio", "Switch_Ratio", "Chimera_Ratio",
        "HPI_Score", "Weighted_Confusion_Index"
    ]
    out_f.write("\t".join(header) + "\n")
    
    global_stats = {"total": 0, "hpi_sum": 0}
    
    for scaff, targets in scaffold_stats.items():
        if not targets: continue
        
        # 3.1 确定 Dominant Identity
        dominant_id = max(targets, key=targets.get)
        dom_chrom, dom_hap = dominant_id
        
        # 3.2 统计 Pure, Switch, Chimera
        l_pure = 0
        l_switch = 0
        l_chimera = 0
        total_len = 0
        
        for (chrom, hap), length in targets.items():
            total_len += length
            if chrom == dom_chrom:
                if hap == dom_hap:
                    l_pure += length
                else:
                    l_switch += length # 同染色体，不同单倍型
            else:
                l_chimera += length # 不同染色体
        
        # 3.3 计算分数
        wci = (args.w_switch * l_switch + args.w_chimera * l_chimera) / total_len
        hpi = max(0.0, 1.0 - wci)
        
        # 记录全局
        global_stats["total"] += total_len
        global_stats["hpi_sum"] += (hpi * total_len)
        
        # 3.4 写入
        dom_str = f"{dom_chrom}|{dom_hap}"
        out_f.write(f"{scaff}\t{dom_str}\t{total_len}\t"
                    f"{l_pure}\t{l_switch}\t{l_chimera}\t"
                    f"{l_pure/total_len:.4f}\t{l_switch/total_len:.4f}\t{l_chimera/total_len:.4f}\t"
                    f"{hpi:.4f}\t{wci:.4f}\n")
                    
    out_f.close()
    
    # 屏幕输出总结
    if global_stats["total"] > 0:
        avg_hpi = global_stats["hpi_sum"] / global_stats["total"]
        print(f"\n[Result] Report saved to: {args.output}", file=sys.stderr)
        print(f"[Result] Global Weighted HPI: {avg_hpi:.4f}", file=sys.stderr)
    else:
        print("\n[Warning] No valid scaffolds evaluated.", file=sys.stderr)

if __name__ == "__main__":
    main()