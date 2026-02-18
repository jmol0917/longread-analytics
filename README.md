# Long-read sequencing analytics portfolio

This repository demonstrates a small, reproducible long-read analytics workflow built as a skills portfolio:
- Read QC
- Alignment to GRCh38 with minimap2
- Structural variant calling with Sniffles2
- Methylation summarisation with modkit
- Reporting (tables + figures)

## Status
What’s implemented so far:
  FASTQ QC + reporting
  FASTQ → BAM alignment pipeline
  SV calling + summary tables + plots
  Modified-base BAM processing + chr1 CpG pileup + QC plots + summary report

## Repo structure
- `data/`
  - `raw/` raw FASTQ/BAM downloads (not tracked in git)
  - `refs/` reference genomes and indexes (not tracked if large)
- `workflow/`
  - `scripts/` executable workflow steps
  - `configs/` config files (paths, threads, references)
  - `logs/` run logs (not tracked)
- `results/`
  - `qc/`, `alignments/`, `sv/`, `methylation/`, `figures/` (not tracked)
- `reports/` short writeups and summaries
- `docs/` tool versions and dataset notes
- `src/` helper Python code for summarisation/plotting
- `notebooks/` exploratory notebooks

## Environment / reproducibility
Create the conda environment:

```bash
conda env create -f environment.yml
conda activate longread-portfolio
```

Tool versions are recorded in docs/versions.txt.

## Workflow overview
1. QC: workflow/scripts/00_qc.sh
2. Align: workflow/scripts/01_align.sh
3. SV calling: workflow/scripts/02_sv_call.sh
4. Methylation: workflow/scripts/03_methylation.sh

## Datasets
Datasets will be documented in docs/datasets.md with accession IDs and download instructions.

## Notes
This repository is intended to demonstrate practical long-read analytics, reproducibility, and clear reporting.
