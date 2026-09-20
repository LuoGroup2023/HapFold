# HapFold

Haplotype-aware assembly and T2T-level scaffolding for diploid genomes.

[![BioConda Downloads](https://img.shields.io/conda/dn/bioconda/hapfold?label=bioconda%20downloads)](https://anaconda.org/bioconda/hapfold)
[![BioConda Version](https://img.shields.io/conda/vn/bioconda/hapfold?label=bioconda)](https://anaconda.org/bioconda/hapfold)
[![License](https://img.shields.io/github/license/LuoGroup2023/HapFold)](https://github.com/LuoGroup2023/HapFold/blob/main/LICENSE)
[![Release](https://img.shields.io/github/v/release/LuoGroup2023/HapFold?sort=semver)](https://github.com/LuoGroup2023/HapFold/releases)
[![Views](https://komarev.com/ghpvc/?username=LuoGroup2023-HapFold&label=Views&color=009E73&style=flat)](https://github.com/LuoGroup2023/HapFold)
[![Stars](https://img.shields.io/github/stars/LuoGroup2023/HapFold?style=social)](https://github.com/LuoGroup2023/HapFold/stargazers)

HapFold combines assembly-graph topology with paired Hi-C or converted Pore-C
contacts. It can run the embedded Hifiasm assembly frontend and HapFold
scaffolding as one workflow, or run mapping and scaffolding separately from
existing Hifiasm GFA files.

## Quick installation

Requirements: Linux, `g++` 8.5 or newer, `make`, `zlib`, and pthreads.

```bash
git clone https://github.com/LuoGroup2023/HapFold.git
cd HapFold
make -j8
```

## Quick start

The shortest complete HiFi + Hi-C command is:

```bash
HapFold run \
  -t THREADS -n CHR_N \
  -1 hic.R1.fastq.gz -2 hic.R2.fastq.gz \
  -o result/asm \
  -- hifi.fastq.gz
```

For ONT simplex reads, pass the native Hifiasm `--ont` option after the
separator:

```bash
HapFold run \
  -t THREADS -n CHR_N \
  -1 hic.R1.fastq.gz -2 hic.R2.fastq.gz \
  -o result/ont \
  -- --ont ont.simplex.fastq.gz
```

`-t` and `-o` are shared workflow options and must occur before `--`.
Everything after `--` is passed to the embedded Hifiasm frontend.
`THREADS` is the number of CPU threads, and `CHR_N` is the expected total
number of chromosomes in the diploid genome (for example, `46` for human).

## Why HapFold?

HapFold is designed to push diploid assemblies from contig-level or partially
scaffolded results toward telomere-to-telomere (T2T)-level reconstruction. Its
main advantages are:

1. **Higher chromosome completeness.** Graph-aware local extension and
   long-range contact scaffolding recover chromosome arms that remain split in
   the initial assembly.
2. **Explicit homolog pairing and phasing.** Bubble relationships, graph
   components and paired-contact evidence are used together so that paternal
   and maternal chromosome paths remain distinguishable during global joining.
3. **More T2T-level assembly products.** HapFold can connect complementary
   T2T-level haplotype paths and can optionally attach real telomeric unitigs
   already present in the assembly graph. It never synthesizes telomere
   sequence.

HapFold is reference-free during assembly and scaffolding. A reference genome
is not required by `run`, `mapping`, or `scaffolding`.

## Commands

```text
HapFold run            embedded Hifiasm followed by HapFold scaffolding
HapFold hifiasm        embedded Hifiasm frontend only
HapFold mapping        contact-read mapping only
HapFold scaffolding    phasing and scaffolding from existing assembly graphs
```

Use `HapFold <command> -h` for the options compiled into the installed binary.

## Complete `run` workflow

```bash
HapFold run \
  -t THREADS \
  -n CHR_N \
  -1 CONTACT_R1.fastq.gz \
  -2 CONTACT_R2.fastq.gz \
  -o OUTPUT_PREFIX \
  [run options] \
  -- [Hifiasm options] ASSEMBLY_READS...
```

In the default `--hifiasm-mode hic` mode, HapFold uses Hifiasm's
exact-extension and chain mapping, supplemented by bounded independent
unique-anchor support. The accepted paired hits are converted directly into
the HapFold contact matrix for downstream phasing and scaffolding.

### `run` options

| Option | Meaning | Default |
|---|---|---|
| `-1`, `--hic1 FILE` | Contact read 1 | required |
| `-2`, `--hic2 FILE` | Contact read 2 | required |
| `-n`, `--chromosomes INT` | Expected diploid chromosome count | required |
| `-t`, `--threads INT` | Threads shared by assembly and scaffolding | `32` |
| `-o`, `--output-prefix STR` | Shared output prefix | `hifiasm.asm` |
| `--hapfold-output DIR` | Final HapFold directory | `<prefix>.hapfold` |
| `--hifiasm-mode hic\|default\|trio` | Embedded assembly mode | `hic` |
| `--high-quality-utg` | Use five correction rounds and enhanced graph cleaning | off |
| `--keep-hifiasm-output` | Retain all raw Hifiasm outputs and restart caches | off |
| `--mapping-output FILE` | Preserved sparse contact matrix | `<output>/mapping.txt` |
| `--utg-gfa FILE` | Override detected `p_utg.gfa` | automatic |
| `-u FILE` | UTG-to-CTG mapping output path | output directory |
| `-d`, `--debug` | Enable debugging output | off |
| `--hic_scaffold_threshold_ratio FLOAT` | Local contact extension threshold | `0.60` |
| `--chain_len_thresh INT` | Long-chain threshold | `12000000` |
| `--scaffold_len_thresh INT` | Direct scaffold-output threshold | `300000` |
| `--global-scaffolding mcl\|legacy` | Global grouping/joining method | `mcl` |
| `--paired-global-merge off\|supported\|inferred` | Homolog-paired global joining | `supported` |
| `--mcl-inflation FLOAT` | MCL inflation; `0` performs an automatic scan | `0` |
| `--forced-chain-links FILE` | Experimental/debug supported-link whitelist | off |
| `--telo` | Enable real telomeric-unitig completion | off |
| `--telo-motif STR` | Telomere motif | `CCCTAA` |

For the complete and current option list, run:

```bash
HapFold run -h -- dummy
HapFold mapping -h
HapFold scaffolding -h
HapFold hifiasm -h
```

### Hybrid Hifiasm mapping options

These native Hifiasm options belong after `--`:

| Option | Meaning | Default |
|---|---|---|
| `--hybrid-min-unique-anchors INT` | Minimum independent unique anchors for a candidate | `1` |
| `--hybrid-unique-weight FLOAT` | Unique-anchor support weight | `0.50` |
| `--hybrid-unique-bonus-cap FLOAT` | Maximum bonus relative to original chain score | `0.25` |

Example with explicit hybrid settings and telomere completion:

```bash
HapFold run \
  -t THREADS -n CHR_N \
  -1 hic.R1.fastq.gz -2 hic.R2.fastq.gz \
  -o result/asm \
  --high-quality-utg \
  --global-scaffolding mcl \
  --paired-global-merge supported \
  --telo --telo-motif CCCTAA \
  -- \
  --dual-scaf --telo-m CCCTAA \
  --hybrid-min-unique-anchors 2 \
  --hybrid-unique-weight 0.5 \
  --hybrid-unique-bonus-cap 0.25 \
  hifi.fastq.gz
```

## Running mapping and scaffolding separately

Extract unitig sequences:

```bash
awk '$1=="S" {print ">"$2"\n"$3}' asm.hic.p_utg.gfa > p_utg.fa
```

Map paired contacts:

```bash
HapFold mapping \
  -t THREADS -k 31 -p 22 \
  -1 hic.R1.fastq.gz -2 hic.R2.fastq.gz \
  -o mapping.txt p_utg.fa
```

For a human-sized graph, `-b 37` enables a `2^37`-bit Bloom filter and can
reduce memory used by the exact k-mer table.

Run phasing and scaffolding:

```bash
HapFold scaffolding \
  -t THREADS -n CHR_N \
  -1 asm.hic.hap1.p_ctg.gfa \
  -2 asm.hic.hap2.p_ctg.gfa \
  --global-scaffolding mcl \
  --paired-global-merge supported \
  mapping.txt asm.hic.p_utg.gfa hapfold_out
```

### `mapping` options

| Option | Meaning | Default |
|---|---|---|
| `-1 FILE`, `-2 FILE` | Paired Hi-C/Pore-C reads | required |
| `-k INT` | k-mer length | `31` |
| `-p INT` | Prefix length | `22` |
| `-b INT` | Bloom-filter size as `2^INT` bits; `0` disables | `0` |
| `-t INT` | Worker threads | `32` |
| `-K INT` | Input chunk size | `100m` |
| `-c FILE` | Optional node-type CSV | none |
| `-L INT` | Mapping/phasing model selector | implementation-specific |
| `-o FILE` | Sparse mapping output | required |

The scaffolding options are the corresponding phasing, global-scaffolding,
paired-merge, telomere and debug options listed for `run` above. In standalone
mode, `-1` and `-2` mean the two Hifiasm haplotype p_ctg GFA files rather than
contact FASTQ files. Standalone scaffolding additionally accepts
`--split-chain-list FILE` to restore named first-contig or internal chain IDs
as source contigs during controlled reruns.

## Pore-C conversion

HapFold accepts Pore-C after multi-contact molecules have been converted into
paired contacts. For a molecule split into fragments A, B and C, the converter
emits AB, AC and BC pairs.

One reproducible route is to generate a monomer BAM with EPI2ME
[`wf-pore-c`](https://github.com/epi2me-labs/wf-pore-c), name-sort it, and run
the bundled converter:

```bash
nextflow run epi2me-labs/wf-pore-c \
  --fastq porec.fastq.gz \
  --ref reference.fa \
  --hi_c true \
  --out_dir wf_pore_c_out

samtools sort -n -@ THREADS \
  -o porec.monomers.name_sorted.bam \
  wf_pore_c_out/<monomer_bam>

g++ -O3 scripts/porec2hic_bam_fast.cpp \
  -o porec2hic_bam_fast -lz

./porec2hic_bam_fast \
  porec.monomers.name_sorted.bam \
  porec.R1.fastq.gz porec.R2.fastq.gz
```

The output pair files can be supplied directly to `-1` and `-2`. Check that
the R1 and R2 record counts are identical before assembly.

## Output files

The principal final files are:

| File | Description |
|---|---|
| `scaffold.fa` | Phased T2T-level scaffolds |
| `hap_contig.fa` | Unresolved/single-chain haplotype contigs |
| `contig_hap_nodes.txt` | Local haplotype membership of graph nodes |
| `contig_composition.txt` | Unitig composition of output contigs/scaffolds |
| `final_sequence_name_map.tsv` | Stable output name to internal path mapping |
| `chromosome_clusters.tsv` | MCL chromosome-group assignments |
| `scaffold_mcl_paths.tsv` | Global scaffold paths selected from clusters |
| `paired_merge_decisions.tsv` | Accepted and rejected paired global joins |
| `telomere_unitig_completion.tsv` | Telomeric-unitig decisions when `--telo` is enabled |
| `phasing_n_gap_restoration.tsv` | Number of unique 100-bp N gaps removed when writing `phasing_hap1.fa` and `phasing_hap2.fa` |
| `n_gap_restoration.tsv` | Per-sequence count and removal audit for unique 100-bp N gaps |

The complete assembly is the union of `scaffold.fa` and `hap_contig.fa`:

```bash
cat hapfold_out/scaffold.fa hapfold_out/hap_contig.fa > hapfold_out/all_scaffold.fa
```

HapFold applies the same conservative restoration rule twice: first when
writing `phasing_hap1.fa` and `phasing_hap2.fa`, and again when writing the
final `scaffold.fa` (and retained `hap_contig.fa`). A 100-bp N gap is removed
only when it is the single N-run and the only ambiguous sequence in that
record. Sequences with multiple gaps, non-100-bp gaps, longer N-runs or a
terminal N-run are left unchanged. Stage counts are recorded in
`phasing_n_gap_restoration.tsv`, and final per-sequence decisions are recorded
in `n_gap_restoration.tsv`.

## Citation

```bibtex
@article{liu2026efficient,
  title   = {Efficient and accurate near telomere-to-telomere haplotype reconstruction of diploid genomes},
  author  = {Liu, Yuansheng and Li, Yichen and Xu, Jialu and Tan, Zhongzheng and Zhang, Wenhai and Wang, Long and Xu, Luohao and Zeng, Xiangxiang and Schoenhuth, Alexander and Luo, Xiao},
  journal = {bioRxiv},
  year    = {2026},
  publisher = {Cold Spring Harbor Laboratory}
}
```

## License

HapFold is distributed under the GNU General Public License v3.0. See
[`LICENSE`](LICENSE).
