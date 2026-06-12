import pandas as pd
import argparse
import gzip

def parse_gff(gff_file):
    """Parse GFF/GTF file and extract coordinate information for all genes"""
    print(f"[Info] Parsing GFF file: {gff_file}")
    genes = []
    open_func = gzip.open if gff_file.endswith('.gz') else open
    
    with open_func(gff_file, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            # Extract information only at the gene level
            if parts[2] != 'gene': continue
            
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]
            
            # Extract Gene Name or Gene ID
            gene_name = "Unknown"
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('Name='):
                    gene_name = attr.split('=')[1]
                elif attr.startswith('gene_name='): # Compatible with GTF format
                    gene_name = attr.split('=')[1].replace('"', '')
                elif attr.startswith('ID=') and gene_name == "Unknown":
                    gene_name = attr.split('=')[1]
                    
            genes.append({
                'Chrom': chrom, 'Gene_Start': start, 'Gene_End': end,
                'Strand': strand, 'Gene_Name': gene_name
            })
    return pd.DataFrame(genes)

def load_detail_tsv(tsv_file):
    """Load the previously generated detail.tsv file"""
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        return df
    except Exception as e:
        print(f"[Error] Failed to read file {tsv_file}: {e}")
        return pd.DataFrame()

def is_gene_fully_covered(chrom, g_start, g_end, detail_df):
    """Check if the gene is continuously and fully covered by a single block of the tool"""
    matches = detail_df[
        (detail_df['T_Name'] == chrom) & 
        (detail_df['T_Start'] <= g_start) & 
        (detail_df['T_End'] >= g_end)
    ]
    return len(matches) > 0

def has_break_near_gene(chrom, g_start, g_end, detail_df, window=500):
    """Check if an assembly break point (T_Start or T_End) appears near the gene"""
    search_start = g_start - window
    search_end = g_end + window
    
    # Search for T_Start or T_End falling within this range
    cond_start = (detail_df['T_Name'] == chrom) & (detail_df['T_Start'] >= search_start) & (detail_df['T_Start'] <= search_end)
    cond_end = (detail_df['T_Name'] == chrom) & (detail_df['T_End'] >= search_start) & (detail_df['T_End'] <= search_end)
    
    return (cond_start | cond_end).any()

def find_rescued_genes(gff_df, baseline_df, target_df, window=500):
    """Core comparison logic: Find genes broken in baseline tool but connected by target tool"""
    rescued_genes = []
    
    total_genes = len(gff_df)
    print(f"[Info] Gene database loaded. Total {total_genes} genes. Starting cross-comparison...")
    
    for idx, row in gff_df.iterrows():
        if idx > 0 and idx % 5000 == 0:
            print(f"  ... Processed {idx}/{total_genes} genes")
            
        chrom = row['Chrom']
        g_start = row['Gene_Start']
        g_end = row['Gene_End']
        
        # 1. Skip if baseline tool already covers it continuously
        if is_gene_fully_covered(chrom, g_start, g_end, baseline_df):
            continue
            
        # 2. Check if baseline tool had a break point near the gene
        if not has_break_near_gene(chrom, g_start, g_end, baseline_df, window):
            continue
            
        # 3. Check if target tool covers it continuously
        if is_gene_fully_covered(chrom, g_start, g_end, target_df):
            rescued_genes.append({
                'Chrom': chrom,
                'Gene_Name': row['Gene_Name'],
                'Gene_Start': g_start,
                'Gene_End': g_end,
                'Strand': row['Strand'],
                'Gene_Length': g_end - g_start
            })
            
    return pd.DataFrame(rescued_genes)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find genes broken in one tool but rescued by another")
    parser.add_argument('-g', '--gff', required=True, help="GFF annotation file")
    parser.add_argument('-b', '--baseline_detail', required=True, help="detail.tsv from baseline tool")
    parser.add_argument('-t', '--target_detail', required=True, help="detail.tsv from target tool")
    parser.add_argument('-o', '--output', default="rescued_genes.tsv", help="Output file path")
    parser.add_argument('-w', '--window', type=int, default=500, help="Window size for breaks (bp)")
    
    args = parser.parse_args()
    
    gff_df = parse_gff(args.gff)
    baseline_df = load_detail_tsv(args.baseline_detail)
    target_df = load_detail_tsv(args.target_detail)
    
    result_df = find_rescued_genes(gff_df, baseline_df, target_df, window=args.window)
    
    if not result_df.empty:
        result_df = result_df.sort_values(by=['Chrom', 'Gene_Start'])
        result_df.to_csv(args.output, sep='\t', index=False)
        print(f"\n[Success] Found {len(result_df)} genes rescued by target tool.")
        print(f"[Output] Results saved to: {args.output}")
    else:
        print("\n[Notice] No matching genes found.")