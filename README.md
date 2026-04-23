# HapFold: Efficient and Accurate Chromosome-Scale Haplotype Reconstruction

[![BioConda Downloads](https://img.shields.io/conda/dn/bioconda/hapfold?label=bioconda%20downloads)](https://anaconda.org/bioconda/hapfold)
[![BioConda Version](https://img.shields.io/conda/vn/bioconda/hapfold?label=bioconda)](https://anaconda.org/bioconda/hapfold)
[![License](https://img.shields.io/github/license/LuoGroup2023/HapFold)](https://github.com/LuoGroup2023/HapFold/blob/main/LICENSE)
[![Release](https://img.shields.io/github/v/release/LuoGroup2023/HapFold?sort=semver)](https://github.com/LuoGroup2023/HapFold/releases)
[![Downloads](https://img.shields.io/github/downloads/LuoGroup2023/HapFold/total)](https://github.com/LuoGroup2023/HapFold/releases)
[![Stars](https://img.shields.io/github/stars/LuoGroup2023/HapFold?style=social)](https://github.com/LuoGroup2023/HapFold/stargazers)

## Description

**HapFold** is a scaffolding framework designed for the highly accurate, chromosome-scale haplotype reconstruction of diploid genomes. 

By uniquely integrating the synergistic features of both **graph-based** and **sequence-based** paradigms, HapFold achieves significantly lower misassignment rates and higher computational efficiency than existing methods. Furthermore, when applied to diploid genomes sequenced with standard ONT simplex reads, HapFold enables the robust and scalable reconstruction of a greater number of **near-T2T (telomere-to-telomere)** assemblies.

## Installation and Dependencies

### Prerequisites
* `g++` (supporting C++9.4 or later)
* `zlib`

### 1. Install via Bioconda (Recommended)
HapFold is officially available on Bioconda. This is the fastest and easiest way to install the tool:

```bash
# Create and activate a new environment
conda create -n hapfold
conda activate hapfold

# Install HapFold
conda install -c bioconda hapfold

# Alternatively, you can use mamba for faster dependency resolution:
# mamba install -c bioconda hapfold
```

### 2. Install from Source Code
If you prefer to compile from source:

```bash
git clone https://github.com/LuoGroup2023/HapFold.git
cd HapFold

# Compile the source code
make -j8
```

## 🚀 Quick Start & Workflow

HapFold utilizes a two-step workflow: Mapping and Resolving. 

### 🎯 Primary Application
HapFold is primarily designed for **diploid genome scaffolding** using Hi-C data. It seamlessly integrates with `hifiasm` outputs and requires three primary GFA files from the initial assembly:
1. Unphased unitig graph (`*.p_utg.gfa`)
2. Haplotype 1 contig graph (`*.hap1.p_ctg.gfa`)
3. Haplotype 2 contig graph (`*.hap2.p_ctg.gfa`)

### General Usage
```text
Usage: HapFold <command> <arguments> <inputs>

Commands:
  resolve_haplotypes    Use Hi-C data to resolve haplotypes and scaffold
  hic_mapping           Map Hi-C data to sequences in the graph
  count                 Count k-mers
  version               Print version number
```

### Step 1: Hi-C Mapping (`hic_mapping`)
Before mapping, you need to extract the node sequences from your hifiasm unitig graph into a FASTA file:
```bash
awk '/^S/{print ">"$2;print $3}' hifiasm_p_utg.gfa > hifiasm_p_utg.fa
```

Then, map the raw Hi-C reads to these node sequences:
```bash
HapFold hic_mapping -t 32 -o map.out hifiasm_p_utg.fa hic.R1.fastq.gz hic.R2.fastq.gz
```

**Key Options for `hic_mapping`:**
* `-t INT`: Number of worker threads [32]
* `-o FILE`: Output file to save the mapping relationships (e.g., `map.out`)
* `-k INT`: k-mer size [31]
* `-p INT`: prefix length [22]
* `-b INT`: set Bloom filter size to `2**INT` bits; `-b 37` is recommended for large/human genomes.

---

### Step 2: Haplotype Resolution (`resolve_haplotypes`)
Once the mapping is complete, use the mapping results alongside the GFA files to resolve haplotypes and build chromosome-scale scaffolds.

```bash
HapFold resolve_haplotypes -t 32 -n chr -u utg_ctg_mappings.csv -i true map.out hifiasm_p_utg.gfa output_dir -1 hap1.p_ctg.gfa -2 hap2.p_ctg.gfa
```
*(Positional arguments: `<mapping_result> <unitig.gfa> <output_directory>`)*

**Key Options for `resolve_haplotypes`:**

| Option | Description |
| :--- | :--- |
| `-t INT` | Number of threads [8]. |
| `-n INT` | Expected number of chromosomes (e.g., `46` for human, `78` for chicken) [0]. |
| `-1 FILE` | **(Required)** Path to haplotype 1 GFA file (`*.hap1.p_ctg.gfa`). |
| `-2 FILE` | **(Required)** Path to haplotype 2 GFA file (`*.hap2.p_ctg.gfa`). |
| `-u FILE` | Path to `utg_to_ctg` relationship file. Highly recommended for accurate graph traversing. |
| `-i BOOL` | Enable identity check on contigs (`true`/`false`) [false]. |
| `-f FILE` | Precomputed identity file path; if omitted but `-i true`, the check will run automatically. |
| `-e STR` | Restriction enzymes separated by comma (e.g., `GATC,GANTC`) [ ]. |
| `-c FILE` | Path to `contig_hap_nodes.txt` (required for specific Hi-C phasing modes). |
| `-p` | Enable plant mode (uses alternative phasing algorithms optimized for complex genomes). |
| `-d`, `--debug` | Enable debug mode to run internal test code functions. |
| `--hic_scaffold_threshold_ratio FLOAT` | Threshold ratio for sequence-based Hi-C scaffolding extensions [0.60]. |
```