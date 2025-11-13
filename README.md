# Base Snakemake for ancient data (eager 3.0.x final_bams/raw/data for all non-UDG treatment and final_bams/trimmed/data)
### Adjustments for runtime:
1) Fill in samples.csv
2) Adjust vcf heterozygosity cutoff (if you want to include lower minor allele freqs for heterozygosity checks), default 0.1
