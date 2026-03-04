# Graph-based Chromosome-Scale Phasing and Scaffolding of Diploid Genomes without Reference Guidance

**GraPhaser** is a novel graph-based algorithm designed to integrate **PacBio HiFi** and **Hi-C** data to produce high-resolution, fully phased haplotypes at the chromosome level. By leveraging the connectivity of assembly graphs, GraPhaser achieves base-level resolution for diploid genomes without the need for a reference genome.

When benchmarking on healthy human genomes (HG002), GraPhaser produced high-quality assemblies with:

* **Continuity:** NG50 > 130 Mb.
* **Accuracy:** Switch/Hamming error rates < 1.5%.
* **Completeness:** > 6.0 Gb.
* **Efficiency:** Complete processing in under 12 hours (an order of magnitude faster than traditional methods).

---

## 🛠 Installation

### Prerequisites

* G++ (supporting C++11 or later)
* Zlib

### Build from Source

```sh
git clone https://github.com/LuoGroup2023/GraPhaser.git
cd GraPhaser && make

```

After compilation, the executable binary `GraPhaser` will be available in the root directory.

---

## 🚀 Execution

GraPhaser follows a two-step workflow to generate fully phased sequences.

### Step 1: Hi-C Mapping to Assembly Graph

First, map the Hi-C reads to the node sequences of the assembly graph (e.g., from `hifiasm`).

> **Note:** You can extract node sequences from a GFA file using:
> `awk '/^S/{print ">"$2;print $3}' hifiasm_r_utg.gfa > hifiasm_r_utg.fa`

```sh
# Map Hi-C reads to graph nodes
GraPhaser hic_mapping -t 32 -o <map.out> <hifiasm_r_utg.fa> <hic.R1.fastq.gz> <hic.R2.fastq.gz>

```

### Step 2: Haplotype Resolution and Scaffolding

Resolve the assembly graph into distinct haplotypes using the mapping results.

```sh
# Resolve haplotypes and produce phased sequences
GraPhaser resolve_haplotypes -t 32 -i true <map.out> <hifiasm_r_utg.gfa> <output_dir>

```

* **Inputs:** `map.out` (from Step 1), `hifiasm_r_utg.gfa` (hifiasm unitig graph).
* **Outputs:** Fully phased sequences in `pred_hap1.fa` and `pred_hap2.fa`.

---

## 📊 Benchmarking Results

The following table demonstrates GraPhaser's performance on the **HG002** dataset using OmniC/Arima genomics data:

| Dataset | Total Size | CPU Time | NG50 | Quality (QV) |
| --- | --- | --- | --- | --- |
| **HG002 (Human)** | ~6.1 Gb | ~5.0 h | **~132 Mb** | **~Q50** |

> For further experiments and reproduction, please refer to `experiments.sh`.

---

## ⚠️ Current Limitations

1. **Centromeric Regions:** Current phased sequences do not yet include centromeric regions.
2. **UL Data:** Integration of Ultra-Long (UL) Nanopore data is not yet implemented.
3. **Trio-mode:** Future versions will include support for trio-hifiasm graphs.

---

<!-- ## 📝 Maintenance

**GraPhaser** is under active development by the **LuoGroup**. We welcome any issues, suggestions, or pull requests.

* **GitHub:** [LuoGroup2023/GraPhaser](https://www.google.com/search?q=https://github.com/LuoGroup2023/GraPhaser)
* **Author:** Yichen-HNU (lyc) -->
