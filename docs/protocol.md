## STEP 1: Setup

1) Established your working environment
	- Confirmed you’re running Ubuntu under WSL2 (Linux environment).
	- Chose a fixed project location:
		/home/jp_molina/projects/longread-analytics

2) Created a clean repo structure
	 - Created a project folder with a standard structure, including:
		- data/ (raw, intermediate, references)
		- results/ (outputs)
		- workflow/ (scripts/configs)
		- docs/ (documentation)

3) Created a README
README describes:
	- what the portfolio is
	- what tools are used
	- what steps the workflow performs
	- how to reproduce it

4) Built a conda environment
	- Name: longread 
	- Installed the core toolchain:
		Python 3.11
		minimap2
		samtools
		sniffles2
		NanoPlot
		modkit (installed via the nanoporetech::modkit channel)
		plus Python analysis tools (pandas, numpy, matplotlib, jupyterlab)
	This is now a reproducible runtime environment.

## STEP 2: Reference Genome setup and QC
1) Downloaded the human reference genome
	- Downloaded the GENCODE GRCh38 primary assembly FASTA into:
		data/refs/

2) Indexed the reference properly
	- Generated two critical index files:
		- minimap2 index (GRCh38.mmi)
		- samtools FASTA index (.fai)
	- This ensures:
		- fast alignment
		- fast downstream access for tools

3) Dataset selection and download
	- Identified a usable ONT human dataset in ONT Open Data
		- Explored the ONT Open Data S3 bucket and discovered:
			- GIAB 2025.01 HG001 data exists, but is mainly stored as:
				- huge aligned BAM files (115–133 GB)
				- raw signal POD5 files (20–35 GB each)
		- Not feasible due to size.
	- Selected a basecalled FASTQ dataset instead (HG002)
		- Located a basecalled “pass” FASTQ dataset:
			- giab_2023.05/flowcells/hg002/.../basecalled/pass/
		- Downloaded all available FASTQ chunks (6 files) into:
			- data/raw/HG002_giab_2023_05/
	- Merged and compressed the dataset for workflow simplicity:
		- HG002_subset.fastq.gz
4) Quality control (QC)
	- Ran NanoPlot QC successfully
		- Commands:
			NanoPlot --fastq data/raw/HG002_giab_2023_05/HG002_subset.fastq.gz -o results/qc/nanoplot_hg002_subset
	- Generated a full QC report bundle
		- The output folder contains:
			- NanoPlot-report.html (main QC dashboard)
			- NanoStats.txt (key summary statistics)
			- read length histograms
			- read length vs quality plots
			- yield plots
	- Interpreted QC:
		- read lengths are strong (N50 ~29 kb)
		- quality is in typical ONT range (mean Q ~12–13)
		- dataset is small (~0.18× human genome coverage), which will mainly affect SV yield
		
## STEP 3: Alignment protocol

1) Input
	- FASTQ: data/raw/HG002_giab_2023_05/HG002_subset.fastq.gz
	- Reference: data/refs/GRCh38.primary_assembly.genome.fa
	- minimap2 index: data/refs/GRCh38.mmi

2) Commands
	- minimap2 -ax map-ont -t 8 data/refs/GRCh38.mmi data/raw/HG002_giab_2023_05/HG002_subset.fastq.gz | samtools sort -@ 8 -o results/			 alignments/HG002_subset.sorted.bam

	- samtools index results/alignments/HG002_subset.sorted.bam

	- samtools flagstat results/alignments/HG002_subset.sorted.bam > results/alignments/HG002_subset.flagstat.txt
