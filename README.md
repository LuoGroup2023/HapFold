# HapFold: Efficient and Accurate T2T-level Haplotype Reconstruction

[![BioConda Downloads](https://img.shields.io/conda/dn/bioconda/hapfold?label=bioconda%20downloads)](https://anaconda.org/bioconda/hapfold)
[![BioConda Version](https://img.shields.io/conda/vn/bioconda/hapfold?label=bioconda)](https://anaconda.org/bioconda/hapfold)
[![License](https://img.shields.io/github/license/LuoGroup2023/HapFold)](https://github.com/LuoGroup2023/HapFold/blob/main/LICENSE)
[![Release](https://img.shields.io/github/v/release/LuoGroup2023/HapFold?sort=semver)](https://github.com/LuoGroup2023/HapFold/releases)
[![Views](https://komarev.com/ghpvc/?username=LuoGroup2023-HapFold&label=Views&color=009E73&style=flat)](https://github.com/LuoGroup2023/HapFold)
[![Stars](https://img.shields.io/github/stars/LuoGroup2023/HapFold?style=social)](https://github.com/LuoGroup2023/HapFold/stargazers)

## Description

**HapFold**  is the first hybrid scaffolding framework that synergistically leverages the complementary strengths of graph-based and sequence-based approaches to achieve chromosome-scale, near-T2T haplotype reconstructions for diploid genomes.

By integrating the topological accuracy of assembly graphs with proximity-guided sequence contiguity, HapFold eliminates the structural errors and chromosomal misassignments common in traditional pipelines. Notably, HapFold accelerates computation by an order of magnitude while delivering superior assembly quality. Even when working with standard Oxford Nanopore Technologies (ONT) simplex reads, it enables the reconstruction of a greater number of near-T2T assemblies, providing a robust and scalable solution for high-fidelity diploid genome assembly.

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
  scaffolding            use Hi-C/Pore-C data to resolve haplotypes
  mapping                map Hi-C/Pore-C data to sequences in the graph
  version                print version number
```

### Step 1: Hi-C Mapping (`hic_mapping`)
Before mapping, you need to extract the node sequences from your hifiasm unitig graph into a FASTA file:
```bash
awk '/^S/{print ">"$2;print $3}' hifiasm_p_utg.gfa > hifiasm_p_utg.fa
```

Then, map the raw Hi-C reads to these node sequences:
```bash
HapFold mapping -t 32 -1 hic.R1.fastq.gz -2 hic.R2.fastq.gz -o mapping.txt hifiasm_p_utg.fa
```

**Key Options for `mapping`:**
* `-1 FILE, -2 FILE `: (Required) Paths to Hi-C forward (R1) and reverse (R2) reads.
* `-t INT`: Number of worker threads [32]
* `-o FILE`: Output file to save the mapping relationships (e.g., `map.out`)
* `-k INT`: k-mer size [31]


---

### Step 2: Haplotype Resolution (`scaffolding`)
Once the mapping is complete, use the mapping results alongside the GFA files to resolve haplotypes and build chromosome-scale scaffolds.

```bash
Usage: HapFold scaffolding [options] <mapping.txt> <assembly.gfa> <output_dir> -1 *.hap1.p_ctg.gfa -2 *.hap2.p_ctg.gfa
```
*(Positional arguments: `<mapping_result> <unitig.gfa> <output_directory>`)*

**Key Options for `scaffolding`:**

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


## Citation
If you use HapFold in your research, please cite:
```text
https://doi.org/10.64898/2026.05.20.726711
```