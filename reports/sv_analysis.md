# Structural Variant Analysis (Oxford Nanopore Long Reads) — HG002 Subset

## Objective
To demonstrate an end-to-end long-read structural variant (SV) workflow, from basecalled Oxford Nanopore reads through alignment and SV detection, using a small HG002 benchmark subset aligned to the GRCh38 human reference genome.

---

## Dataset
- Sample: HG002 (Genome in a Bottle benchmark sample)
- Sequencing technology: Oxford Nanopore Technologies (ONT)
- Input format: basecalled FASTQ (pass reads)
- Total reads: 24,000
- Total bases: ~541 Mb
- Approximate genome coverage: ~0.18× (human genome ~3.0 Gb)

---

## Workflow steps

### 1) Quality control
QC was performed using NanoPlot to characterise read length and basecall quality.

Key observations:
- Mean read length: ~22.6 kb
- Read length N50: ~29.0 kb
- Mean read quality: ~Q12.5

### 2) Alignment to GRCh38
Reads were aligned to GRCh38 using minimap2 (ONT preset), followed by BAM sorting, indexing, and alignment summarisation using samtools.

Alignment metrics:
- Primary reads: 24,000
- Primary mapped: 23,941 (~99.75%)
- Supplementary alignments: 3,202
- Secondary alignments: 7,185

Interpretation:
High mapping rate indicates correct reference usage and strong data integrity. Supplementary alignments suggest breakpoint-spanning reads, supporting structural variant detection.

### 3) Structural variant calling
Structural variants were called using Sniffles2 from the coordinate-sorted, indexed BAM file.

Output:
- `results/sv/HG002_subset.sv.vcf` (SV callset)

---

## Key results

### Structural variant type distribution
- Total SV calls: 71
- PASS variants: 67
- Non-PASS variants: 4

SVTYPE distribution:
- Deletions (DEL): 44
- Insertions (INS): 26
- Duplications (DUP): 1
- Inversions (INV): 0
- Breakends / translocations (BND): 0

### Structural variant size distribution (absolute SVLEN)
- Min: 49 bp
- Median: 255 bp
- Mean: 926 bp
- Max: 22,724 bp (~22.7 kb)

Interpretation:
Most variants are in the small-to-medium SV size range (hundreds of bases), with at least one larger event (~22.7 kb). This demonstrates long-read capacity for detecting SVs spanning from tens of bases to tens of kilobases.

---

## Key outputs (files)
- `results/sv/HG002_subset.sv.vcf` — structural variant calls (VCF)
- `results/sv/summary.tsv` — SVTYPE counts + SVLEN summary statistics
- `results/sv/sv_summary.png` — SV length histogram

---

## Tool versions

- minimap2: 2.30-r1287
- samtools: 1.23
- sniffles: 2.7.2
- NanoPlot: 1.46.2
- python: 3.11.14

---

## Limitations (important)
- Coverage is very low (~0.18×), reducing sensitivity for SV detection.
- No benchmarking against GIAB truth sets was performed.
- SV calls in repetitive regions may include false positives and haven't been manually curated yet.
- Inversions and translocations were not detected in this subset, likely due to limited coverage and read support.

---

## Translation to cancer genomics
Long-read SV workflows are directly relevant to cancer genomics because many clinically and biologically important events are structural rather than single-base changes. Long reads improve detection and interpretation of:

- Complex rearrangements and clustered breakpoints
- Large insertions/deletions disrupting tumour suppressors or regulatory elements
- Structural drivers underlying gene fusions
- Rearrangements spanning repetitive or low-complexity regions

In tumour samples, SVs often occur in subclonal populations. Long reads can support tumour heterogeneity analysis by providing breakpoint-level evidence in single molecules, although robust quantification requires higher coverage.


