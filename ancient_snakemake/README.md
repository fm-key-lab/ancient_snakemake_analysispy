# Base Snakemake for ancient data (eager 3.0.x final_bams/raw/data for all non-UDG treatment and final_bams/trimmed/data)
### Adjustments for runtime:
1) Fill in samples.csv
2) Adjust vcf heterozygosity cutoff in posteager.Snakefile (if you want to include potential SNPs at a lower (or higher) minor allele freqs, eg for heterozygosity checks), default 0.1, if not running heterozygosity checks set this to 0.9
3) Run (bash snakemakeslurm.sh --posteager or --cmt for running only the first or second half, if errors in only one)
