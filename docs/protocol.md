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
	
3) Outputs
	- HG002_subset.samtools_stats.txt:
		Command:
			samtools stats results/alignments/HG002_subset.sorted.bam > results/alignments/HG002_subset.samtools_stats.txt

		Alignment metrics:
			- Primary reads: 24,000
			- Primary mapped: 23,941 (~99.75%)
			- Supplementary alignments: 3,202
			- Secondary alignments: 7,185

## STEP 4: Structureal Variant Calling (Sniffles2)

1) Input
	- Aligned BAM: results/alignments/HG002_subset.sorted.bam
	- Sniffles2 also uses: results/alignments/HG002_subset.sorted.bam.bai
			- For access to genomic regions, fast jumping to coordinates, scanning each chromosome
	- Reference genome: data/refs/GRCh38.primary_assembly.genome.fa

2) Commands	
	- Ran Sniffles2 for structural variant calling:
			- sniffles --input results/alignments/HG002_subset.sorted.bam --vcf results/sv/HG002_subset.sv.vcf --threads 8
	- Output: HG002_subset.sv.vcf (Variant Call Format file containing structural variants)
	
3) Outputs
	- Variant summary statistics
		- Counted total structural variants:
			- grep -v "^#" results/sv/HG002_subset.sv.vcf | wc -l
			- Total calls: 71
	
		- Filter:
			- grep -v "^#" results/sv/HG002_subset.sv.vcf | cut -f7 | sort | uniq -c
			- Output:
				- 67 PASS  
				- 4 GT (lower-confidence calls)
				
	- Structural variant type distribution
		- Extracted SVTYPE field:
			grep -v "^#" results/sv/HG002_subset.sv.vcf | cut -f8 | grep -o "SVTYPE=[A-Z]*" | sort | uniq -c | sort -nr
		- Output: Distribution
			44 deletions (DEL)
			26 insertions (INS)
			1 duplication (DUP)
			
	- Structural variant size statistics
		- Total variants: 71
		- Minimum absolute SV length: 49 bp
		- Median absolute SV length: 255 bp
		- Mean absolute SV length: 926 bp
		- Maximum absolute SV length: 22,724 bp
		- Interpretation:
			- Majority of variants are small-to-medium scale (hundreds of bases).
			- Presence of large events (~22.7 kb) demonstrates long-read capability for detecting sizeable genomic rearrangements.
			- Despite low genome coverage (~0.18×), the high mapping rate (~99.7%) enabled robust detection of structural variant signals.
			
	- Reports, tables and figures
		- results/sv/HG002_subset.sv.vcf (SV callset - .vcf stands for Variant Call Format)
		- results/sv/summary.tsv (counts + size distribution)
		- results/sv/sv_summary.png (SV length histogram)
		
## STEP 5: SV reporting

1) Created report file
	Created the report directory:
		- reports/sv_analysis.md
		
	Contents:
	- Workflow steps
		QC (NanoPlot)
		Alignment (minimap2 + samtools)
		Structural variant calling (Sniffles2)
		SV summarisation (summary.tsv + sv_summary.png)

	- Tool versions
		minimap2
		samtools
		sniffles2
		NanoPlot
		python

	- Key outputs
		results/sv/HG002_subset.sv.vcf
		results/sv/summary.tsv
		results/sv/sv_summary.png
	
	- Key results
	
	- Limitations
		Very low coverage (~0.18×), reducing sensitivity
		No benchmarking against GIAB truth set
		No manual curation of calls
		No short-read integration
		SV types like INV/BND not detected in this subset, likely due to limited read support

	- Translation to cancer genomics
		relevance of SVs to cancer drivers
		complex rearrangements
		gene fusions
		tumour heterogeneity
		
## STEP 6: Methylation / Epigenomics

1) Inputs
	ModBAM (contains MM/ML tags):
		data/raw/modbases_big/20210511_1515_X2_FAQ32637_9b683def.bam

	Updated ModBAM (clean MM/ML tags):
		results/modbases/20210511_1515_X2_FAQ32637_9b683def.MMML.bam

	Reference FASTA:
		data/refs/GRCh38.primary_assembly.genome.fa

	FASTA index:
		data/refs/GRCh38.primary_assembly.genome.fa.fai
		
2) Commands
	Download raw BAM
		BAM_S3="s3://ont-open-data/gm24385_mod_2021.09/extra_analysis/alignment/20210511_1515_X2_FAQ32637_9b683def.bam"
		BAM_LOCAL="data/raw/modbases_big/$(basename "$BAM_S3")"
		aws s3 cp --no-sign-request "$BAM_S3" "$BAM_LOCAL"
	
	Update tags from raw BAM
		modkit update-tags "$BAM_LOCAL" results/modbases/$(basename "$BAM_LOCAL" .bam).MMML.bam
	
	Confirm BAM contains modified base tags
		BAM_MM="results/modbases/20210511_1515_X2_FAQ32637_9b683def.MMML.bam"

		samtools view "$BAM_MM" | head -n 2000 | grep -c $'\tMM:Z:'
		samtools view "$BAM_MM" | head -n 2000 | grep -c $'\tML:B:C,'
		
	Generate a chr1-only modkit summary
		modkit summary --region chr1 --threads 1 "$BAM_MM" > results/modbases/20210511_1515_X2_FAQ32637_9b683def.MMML.chr1.modkit_summary.txt
		
		Output produced:
			results/modbases/20210511_1515_X2_FAQ32637_9b683def.MMML.chr1.modkit_summary.txt
			
	Export CpG methylation calls for chr1
		modkit pileup --region chr1 --cpg --ref data/refs/GRCh38.primary_assembly.genome.fa "$BAM_MM" results/modbases/20210511_1515_X2_FAQ32637_9b683def.MMML.chr1.cpg.bed
		
		then compressed it:
			gzip -f results/modbases/20210511_1515_X2_FAQ32637_9b683def.MMML.chr1.cpg.bed
			
	Generate chr1 methylation QC summary files
		- File 1 Convert modkit pileup BED → compact TSV (pos, strand, cov, beta)
			IN_BED_GZ="results/modbases/20210511_1515_X2_FAQ32637_9b683def.MMML.chr1.cpg.bed.gz"
			OUT_TSV="results/modbases/chr1.cpg.beta_cov.tsv.gz"

			zcat "$IN_BED_GZ" awk 'BEGIN{OFS="\t"}{pos=$2+1; cov=$10; beta=$11; print $1,pos,$6,cov,beta}' gzip -c > "$OUT_TSV"
			ls -lh "$OUT_TSV"
		
			Ouptut: results/modbases/chr1.cpg.beta_cov.tsv.gz
		
		- File 2 Plot beta histogram
			OUT_PNG="results/modbases/chr1_beta_hist.png"

			python3 - <<'PY'
			import gzip, numpy as np, matplotlib.pyplot as plt

			tsv = "results/modbases/chr1.cpg.beta_cov.tsv.gz"
			betas = []
			with gzip.open(tsv, "rt") as f:
				for line in f:
					_,_,_,cov,beta = line.rstrip("\n").split("\t")
					# optional: ignore singletons
					# if int(cov) < 2: continue
					b = float(beta)
					if 0.0 <= b <= 1.0:
						betas.append(b)

			betas = np.array(betas, dtype=float)
			plt.figure()
			plt.hist(betas, bins=50)
			plt.xlabel("CpG methylation fraction (beta)")
			plt.ylabel("Count (sites)")
			plt.title("chr1 CpG methylation fraction distribution")
			plt.tight_layout()
			plt.savefig("results/modbases/chr1_beta_hist.png", dpi=200)
			print("Wrote:", "results/modbases/chr1_beta_hist.png", "n=", betas.size)
			PY

			ls -lh "$OUT_PNG"
			
			Output: results/modbases/chr1_beta_hist.png
			
		- File 3 Plot coverage vs beta scatter
			OUT_PNG="results/modbases/chr1_cov_vs_beta.png"

			python3 - <<'PY'
			import gzip, random
			import numpy as np
			import matplotlib.pyplot as plt

			tsv = "results/modbases/chr1.cpg.beta_cov.tsv.gz"
			covs = []
			betas = []

			# downsample to keep plot light
			p_keep = 0.02  # keep ~2% of sites; adjust 0.005–0.05 if desired

			with gzip.open(tsv, "rt") as f:
				for line in f:
					if random.random() > p_keep:
						continue
					_,_,_,cov,beta = line.rstrip("\n").split("\t")
					c = int(cov); b = float(beta)
					if c <= 0:
						continue
					covs.append(c)
					betas.append(b)

			covs = np.array(covs)
			betas = np.array(betas)

			plt.figure()
			plt.scatter(covs, betas, s=2)
			ls -lh "$OUT_PNG"results/modbases/chr1_cov_vs_beta.png", "n=", covs.size)
			
			Output: results/modbases/chr1_cov_vs_beta.png
			
		- File 4 Generate reports/qc_summary_chr1.txt from the TSV
			OUT_TXT="reports/qc_summary_chr1.txt"
			
			python3 - <<'PY'
			import gzip, numpy as np

			tsv = "results/modbases/chr1.cpg.beta_cov.tsv.gz"
			cov=[]; beta=[]
			with gzip.open(tsv,"rt") as f:
				for line in f:
					_,_,_,c,b = line.rstrip("\n").split("\t")
					cov.append(int(c)); beta.append(float(b))
			cov=np.array(cov); beta=np.array(beta)

			def pct(x): return float(np.percentile(x, [1,5,25,50,75,95,99])[3])

			with open("reports/qc_summary_chr1.txt","w") as out:
				out.write(f"Sites (chr1 CpG): {beta.size}\n")
				out.write(f"Coverage median: {np.median(cov):.2f}\n")
				out.write(f"Coverage mean: {np.mean(cov):.2f}\n")
				out.write(f"Beta mean: {np.mean(beta):.4f}\n")
				out.write(f"Beta median: {np.median(beta):.4f}\n")
				out.write("Beta percentiles (1,5,25,50,75,95,99): " +
						  ",".join([f'{p:.4f}' for p in np.percentile(beta,[1,5,25,50,75,95,99])]) + "\n")
			print("Wrote reports/qc_summary_chr1.txt")
			PY

			cat "$OUT_TXT"
			
			Output: reports/qc_summary_chr1.txt
			
3) Outputs and interpretation
	Key metrics:
		Sites (chr1 CpG): 1,397,156
		Coverage median: 1.00
		Coverage mean: 1.60
		Beta mean: 0.0000
		Beta median: 0.0000
		
	Interpretation:
		Coverage is extremely shallow (mostly 1×), so methylation estimates are noisy.
		The near-zero beta values strongly suggest either:
			the modification type is not actually 5mC at CpG, OR
			the tags represent probabilities but thresholding removed calls, OR
			the dataset is not truly “methylation-called” in the way we assumed (even though MM/ML exist)