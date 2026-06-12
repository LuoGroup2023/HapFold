#!/bin/bash
set -euo pipefail

REF="ref.fasta"
ASSEMBLY_FILE="assembly.fa"
PAF_DIR="alignment.paf"
GENOME_INFO="genome_info.txt"
THREADS=16

command -v samtools >/dev/null 2>&1 || { echo "Error: samtools required but not installed."; exit 1; }

echo "Generating $GENOME_INFO from $REF..."
samtools faidx "$REF"
cut -f1,2 "${REF}.fai" | awk '{print $1 " (total length: " $2 " bp)"}' > "$GENOME_INFO"

echo "Running minimap2 with $THREADS threads..."
[[ ! -f "$REF" ]] && { echo "Error: $REF not found."; exit 1; }
[[ ! -f "$ASSEMBLY_FILE" ]] && { echo "Error: $ASSEMBLY_FILE not found."; exit 1; }

minimap2 -x asm5 -t "$THREADS" --secondary=no -N 0 -c "$REF" "$ASSEMBLY_FILE" > "$PAF_DIR"

echo "Running evaluation..."
if [[ -f "eval.py" ]]; then
    python3 eval.py -i "$PAF_DIR" -g "$GENOME_INFO"
else
    echo "Error: eval.py not found in current directory."
    exit 1
fi

echo "All tasks completed successfully."