# Base Snakemake for ancient data (eager 3.0.x final_bams/raw/data for all non-UDG treatment and final_bams/trimmed/data)
### Adjustments for runtime:
1) Fill in samples.csv
2) Make optional adjustments:
   1) Adjust vcf heterozygosity cutoff in posteager.Snakefile (if you want to include potential SNPs at a lower (or higher) minor allele freqs, eg for heterozygosity checks), default 0.1, if not running heterozygosity checks set this to 0.9
   2) Add a bedfile (regions.bed) for creating CMT only across positions within bedfile (eg if only want genic regions, or genes within core genome) -- this can be also do by filtering from called SNP positions at local analysis step, but this can be problematic due to mapping issues
3) Run (bash snakemakeslurm.sh --posteager or --cmt for running only the first or second half, if errors in only one (eg posteager OK but cmt not))
