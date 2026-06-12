import pandas as pd
import argparse
import re
import os

def get_base_chrom(t_name):
    base = re.sub(r'MATERNAL|PATERNAL|MATERAL|PATERAL', '', t_name, flags=re.IGNORECASE)
    base = base.strip('_')  
    return base

def get_haplotype(t_name):
    t_name = t_name.upper()
    if 'MATERNAL' in t_name:
        return 'MAT'
    elif 'PATERNAL' in t_name:
        return 'PAT'
    else:
        return 'UNK'

def parse_genome_info(genome_info_file):
    chrom_lengths = {}
    if not os.path.exists(genome_info_file):
        return chrom_lengths
    with open(genome_info_file, 'r', encoding='utf-8') as f:
        for line in f:
            if '(' in line and 'total length:' in line:
                match = re.search(r'(\S+)\s+\(total length:\s+(\d+)\s+bp', line)
                if match:
                    chrom_lengths[match.group(1)] = int(match.group(2))
    return chrom_lengths

def merge_intervals(intervals):
    if not intervals: return []
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def get_overlap_len(s1, e1, s2, e2):
    return max(0, min(e1, e2) - max(s1, s2))

def calculate_igi_components_row_by_row(blocks, gap_limit=1000000):
    if blocks.empty:
        return 0.0, 0.0, 0.0, 0.0, 0.0
        
    total_len = blocks['Block_Len'].sum()
    if total_len == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0
        
    strand_sums = blocks.groupby('Strand')['Block_Len'].sum()
    primary_strand = strand_sums.idxmax()
    
    covered_refs = []
    prev_valid_t_end = None
    prev_valid_t_start = None
    
    len_dup = 0.0
    len_inv = 0.0
    len_reloc = 0.0
    len_cor = 0.0
    
    blocks = blocks.sort_values('Q_Start')
    
    for _, row in blocks.iterrows():
        b_len = row['Block_Len']
        t_start, t_end = row['T_Start'], row['T_End']
        strand = row['Strand']
        
        overlap_ref = 0
        merged_covered = merge_intervals(covered_refs)
        for cs, ce in merged_covered:
            overlap_ref += get_overlap_len(t_start, t_end, cs, ce)
            
        ref_len = t_end - t_start
        dup_ratio = overlap_ref / ref_len if ref_len > 0 else 0
        dup_ratio = min(1.0, dup_ratio)
        
        dup_L = b_len * dup_ratio
        rem_L = b_len - dup_L
        
        covered_refs.append((t_start, t_end))
        
        if strand != primary_strand:
            inv_L = rem_L
            reloc_L = 0.0
            cor_L = 0.0
        else:
            inv_L = 0.0
            rem_L2 = rem_L
            
            reloc_L = 0.0
            if prev_valid_t_end is not None:
                if primary_strand == '+':
                    dist = t_start - prev_valid_t_end
                else:
                    dist = prev_valid_t_start - t_end
                    
                if dist < -0.1 * ref_len or dist > gap_limit:
                    reloc_L = rem_L2
            
            cor_L = rem_L2 - reloc_L
            
            prev_valid_t_end = t_end
            prev_valid_t_start = t_start
            
        len_dup += dup_L
        len_inv += inv_L
        len_reloc += reloc_L
        len_cor += cor_L
        
    p_dup = len_dup / total_len
    p_inv = len_inv / total_len
    p_reloc = len_reloc / total_len
    p_cor = len_cor / total_len
    
    total_p = p_dup + p_inv + p_reloc + p_cor
    if total_p > 0:
        p_dup /= total_p
        p_inv /= total_p
        p_reloc /= total_p
        p_cor /= total_p
        
    weighted_igi = 0.6 * p_dup + 0.3 * p_inv + 0.1 * p_reloc
    
    return p_inv, p_dup, p_reloc, p_cor, weighted_igi


def run_assembly_evaluator(args):
    chrom_lengths = parse_genome_info(args.genome)
    paf_cols = ['Q_Name', 'Q_Len', 'Q_Start', 'Q_End', 'Strand', 'T_Name', 'T_Len', 'T_Start', 'T_End', 'Matches', 'Block_Len', 'MapQ']
    try:
        df = pd.read_csv(args.input, sep='\t', header=None, usecols=range(12), names=paf_cols)
    except Exception as e:
        print(f"load PAF fail: {e}")
        return

    global_stats = {
        'total_aligned_len': 0, 'total_primary_hap_len': 0,
        'weighted_cgi_sum': 0, 'weighted_hgi_sum': 0, 
        'weighted_igi_sum': 0
    }

    contig_scores_list = []
    mapping_summary_list = []
    all_filtered_dfs = []

    print(f"Running Row-by-Row Assessment for IGI (Duplication, Inversion, Relocation)...")
    
    for q_name, group in df.groupby('Q_Name'):
        q_total_len = group['Q_Len'].iloc[0]
        
        # --- A. Tiling Path ---
        sorted_aligns = group.sort_values(by='Matches', ascending=False)
        selected_tiles, covered_intervals = [], []
        for _, row in sorted_aligns.iterrows():
            qs, qe = row['Q_Start'], row['Q_End']
            overlap = sum(max(0, min(qe, e) - max(qs, s)) for s, e in covered_intervals)
            if overlap < ((qe - qs) * 0.1): 
                selected_tiles.append(row)
                covered_intervals.append((qs, qe))
        
        if not selected_tiles: continue
        
        raw_tiled_df = pd.DataFrame(selected_tiles).sort_values(by='Q_Start')
        tiled_df = raw_tiled_df[raw_tiled_df['Block_Len'] >= args.length].copy()
        if tiled_df.empty: 
            continue
        all_filtered_dfs.append(tiled_df)

        # tiled_df['Ref_Chrom_Base'] = tiled_df['T_Name'].apply(lambda x: x.split('_')[0])
        tiled_df['Ref_Chrom_Base'] = tiled_df['T_Name'].apply(get_base_chrom)
        tiled_df['Haplotype'] = tiled_df['T_Name'].apply(get_haplotype)
        chrom_sums = tiled_df.groupby('Ref_Chrom_Base')['Block_Len'].sum()
        primary_prefix = chrom_sums.idxmax()
        contig_aligned_len = chrom_sums.sum()
        for t_name, mapped_len in tiled_df.groupby('T_Name')['Block_Len'].sum().items():
            ref_len = chrom_lengths.get(t_name, 1)
            frac = mapped_len / ref_len
            if frac >= args.threshold:
                mapping_summary_list.append({
                    'Contig': q_name, 'Ref_Chrom': t_name, 'Mapped_Length': mapped_len,
                    'Genome_Fraction': round(frac, 6), 'Is_Primary': (t_name.split('_')[0] == primary_prefix)
                })

        # [1] CGI
        if contig_aligned_len > 0:
            sum_p2 = sum((length / contig_aligned_len)**2 for length in chrom_sums)
            cgi_score = 1.0 - sum_p2
        else:
            cgi_score = 0.0
        global_stats['total_aligned_len'] += contig_aligned_len
        global_stats['weighted_cgi_sum'] += (cgi_score * contig_aligned_len)

        # [2] HGI
        primary_blocks = tiled_df[tiled_df['Ref_Chrom_Base'] == primary_prefix].copy()
        total_hap_len, hgi_score = 0, 0.0
        if not primary_blocks.empty:
            mat_len = primary_blocks[primary_blocks['T_Name'].str.contains('MAT', case=False)]['Block_Len'].sum()
            pat_len = primary_blocks[primary_blocks['T_Name'].str.contains('PAT', case=False)]['Block_Len'].sum()
            if mat_len >= pat_len:
                primary_map_len = mat_len
                other_map_len = pat_len
            else:
                primary_map_len = pat_len
                other_map_len = mat_len
            total_hap_len = mat_len + pat_len
            if total_hap_len > 0:
                p_mat = mat_len / total_hap_len
                p_pat = pat_len / total_hap_len
                hgi_score = (1.0 - (p_mat**2 + p_pat**2)) * 2.0 
        global_stats['total_primary_hap_len'] += total_hap_len
        global_stats['weighted_hgi_sum'] += (hgi_score * total_hap_len)

        # [3] IGI 
        if not primary_blocks.empty:
            p_inv, p_red, p_reloc, p_cor, igi_score = calculate_igi_components_row_by_row(primary_blocks, gap_limit=1000000)
        else:
            p_inv, p_red, p_reloc, p_cor, igi_score = 0.0, 0.0, 0.0, 0.0, 0.0
            
        global_stats['weighted_igi_sum'] += (igi_score * contig_aligned_len)

        total_error_sum = cgi_score + hgi_score + igi_score
        contig_scores_list.append({
            'Contig': q_name, 'Contig_Len': q_total_len, 'Primary_Ref': primary_prefix,
            'Primary_Map_Len': primary_map_len, 'Other_Hap_Len': other_map_len,
            'Total_Error_Sum': round(total_error_sum, 4), 'CGI': round(cgi_score, 4),
            'HGI': round(hgi_score, 4), 'IGI': round(igi_score, 4), 
            'P_Duplication': round(p_red, 4), 'P_Inversion': round(p_inv, 4), 
            'P_Relocation': round(p_reloc, 4), 'P_Correct': round(p_cor, 4)
        })

    if not contig_scores_list: return
    scores_df = pd.DataFrame(contig_scores_list)
    output_cols = ['Contig', 'Contig_Len', 'Primary_Ref', 'Primary_Map_Len', 'Other_Hap_Len', 'Total_Error_Sum', 'CGI', 'HGI', 'IGI', 
                   'P_Duplication', 'P_Inversion', 'P_Relocation', 'P_Correct']
    scores_df = scores_df[output_cols].sort_values(by='Total_Error_Sum', ascending=False)
    scores_df.to_csv(args.output_contig_error, sep='\t', index=False)
    print(f"[Output] Contig Errors: {args.output_contig_error}")

    if all_filtered_dfs: pd.concat(all_filtered_dfs).to_csv(args.output_detail, sep='\t', index=False)
    if mapping_summary_list: 
        pd.DataFrame(mapping_summary_list).to_csv(args.output_summary, sep='\t', index=False)
        print(f"[Output] Mapping Summary: {args.output_summary}")
    global_assess = []
    total_len = global_stats['total_aligned_len']
    if total_len > 0:
        g_cgi = global_stats['weighted_cgi_sum'] / total_len
        g_igi = global_stats['weighted_igi_sum'] / total_len
    else:
        g_cgi, g_igi = 0.0, 0.0

    total_hap_len = global_stats['total_primary_hap_len']
    g_hgi = (global_stats['weighted_hgi_sum'] / total_hap_len) if total_hap_len > 0 else 0.0

    g_total = g_cgi + g_hgi + g_igi
    global_assess.append({'Metric': 'Global_Total_Error_Sum', 'Value': round(g_total, 6)})
    global_assess.append({'Metric': 'Global_CGI', 'Value': round(g_cgi, 6)})
    global_assess.append({'Metric': 'Global_HGI', 'Value': round(g_hgi, 6)})
    global_assess.append({'Metric': 'Global_IGI', 'Value': round(g_igi, 6)})

    pd.DataFrame(global_assess).to_csv(args.output_assess, sep='\t', index=False)
    print(f"[Output] Global Assessment: {args.output_assess}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Nature-Tier Scaffolding Evaluator - IGI Updated Edition')
    parser.add_argument('-i', '--input', required=True, help='Input PAF file')
    parser.add_argument('-g', '--genome', required=True, help='genome_info.txt')
    parser.add_argument('-l', '--length', type=int, default=500000, help='Length Filter (bp)')
    parser.add_argument('-d', '--output_detail', default='detail.tsv')
    parser.add_argument('-s', '--output_summary', default='mapping_summary.tsv')
    parser.add_argument('-e', '--output_contig_error', default='error.tsv')
    parser.add_argument('-a', '--output_assess', default='assess.tsv')
    parser.add_argument('-t', '--threshold', type=float, default=0.1) 
    args = parser.parse_args()
    run_assembly_evaluator(args)