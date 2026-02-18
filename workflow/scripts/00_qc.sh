#!/usr/bin/env bash
set -euo pipefail

FASTQ="data/raw/HG002_giab_2023_05/HG002_subset.fastq.gz"
OUTDIR="results/qc/nanoplot_hg002_subset"

echo "Running NanoPlot QC..."

mkdir -p "$OUTDIR"

NanoPlot \
  --fastq "$FASTQ" \
  --outdir "$OUTDIR"

echo "QC complete."
echo "Output: $OUTDIR"

