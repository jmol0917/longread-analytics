#!/usr/bin/env bash
set -euo pipefail

FASTQ="data/raw/HG002_giab_2023_05/HG002_subset.fastq.gz"
REF_INDEX="data/refs/GRCh38.mmi"
OUTDIR="results/alignments"
BAM="$OUTDIR/HG002_subset.sorted.bam"

echo "Running alignment..."

mkdir -p "$OUTDIR"

minimap2 -ax map-ont -t 8 "$REF_INDEX" "$FASTQ" \
  | samtools sort -@ 8 -o "$BAM"

samtools index "$BAM"

samtools flagstat "$BAM" > "$OUTDIR/HG002_subset.flagstat.txt"

echo "Alignment complete."
echo "Output BAM: $BAM"

