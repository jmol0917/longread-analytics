#!/usr/bin/env bash
set -euo pipefail

BAM="results/alignments/HG002_subset.sorted.bam"
OUTDIR="results/sv"
VCF="$OUTDIR/HG002_subset.sv.vcf"

echo "Running Sniffles2 SV calling..."

mkdir -p "$OUTDIR"

sniffles \
  --input "$BAM" \
  --vcf "$VCF" \
  --threads 8

echo "SV calling complete."
echo "Output VCF: $VCF"

