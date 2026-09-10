#############################################
# KEYLAB ANCIENT SNAKEFILE FOR MAPPING STEP #
#############################################
''' PRE-SNAKEMAKE '''

import sys
import os
SCRIPTS_DIRECTORY = "./scripts"
sys.path.insert(0, SCRIPTS_DIRECTORY)

# from import read_samplesCSV # not needed since this function is in the next file
from read_move_link_samplesCSV import *

minMAF = 0.1

## Format: Path,Sample,ReferenceGenome,ProviderName,Subject
spls = "samples.csv"
[PATH_ls,SAMPLE_ls,REF_Genome_ls,CALLINDELS_ls,OUTGROUP_ls] = read_samplesCSV(spls)
bams_ls = get_bams(SAMPLE_ls, REF_Genome_ls)
[REF_Genome_ext_ls, SAMPLE_ext_ls] = parse_multi_genome_smpls(SAMPLE_ls, REF_Genome_ls)
sample_to_reference = dict(zip(SAMPLE_ls, REF_Genome_ls))

# grab current working directory for qc rules to use
current_directory = os.getcwd()


''' SNAKEMAKE '''

rule all:
  input:
    expand("2-quals/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.quals.npz", zip, sampleID=SAMPLE_ls, reference=REF_Genome_ls),
    expand("3-diversity/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.diversity.npz", zip, sampleID=SAMPLE_ls, reference=REF_Genome_ls),
    expand("1-vcf/ref_{reference}_freebayes_raw_joint_calls.vcf",reference=set(REF_Genome_ext_ls)),
    expand("4-bed_files/{sampleID}_genome_coverage_hist.tsv.gz", sampleID=SAMPLE_ls),
    expand("4-bed_files/{sampleID}_merged_zero_covg_regions.tsv.gz", sampleID=SAMPLE_ls),
    "samples_case.csv",
    "cleanUp_done.txt",
    "samples.csv"

rule freebayes_indels:
  input:
    non_outgroup_bam_list="0-freebayes_input/ref_{reference}_non_outgroup_bams.txt", 
    fai="/nexus/posix0/MPIIB-keylab/reference_genomes/{reference}/genome.fasta.fai",
    ref="/nexus/posix0/MPIIB-keylab/reference_genomes/{reference}/genome.fasta",
  output:
    vcf_indels="1-vcf/ref_{reference}_non_outgroup_indels_complex.vcf.gz",
    vcf_raw="1-vcf/ref_{reference}_freebayes_raw_joint_calls.vcf",
  params:
    regions = "regions.bed",
  conda:
    "envs/freebayes.yaml", 
  shell:
    """
        if [ ! -s {input.non_outgroup_bam_list} ]; then
            > {output.vcf_raw} ;
            > {output.vcf_indels} ;
        elif [ -e {params.regions} ]; then 
            freebayes 72 -t {params.regions} -f {input.ref} -p 1 -L {input.non_outgroup_bam_list} > {output.vcf_raw} ;
            egrep '#|ins|del|complex' {output.vcf_raw} | gzip -c > {output.vcf_indels} ;
        else
            freebayes-parallel <(fasta_generate_regions.py {input.fai} 100000) 72 -f {input.ref} -p 1 -L {input.non_outgroup_bam_list} > {output.vcf_raw} ;
            egrep '#|ins|del|complex' {output.vcf_raw} | gzip -c > {output.vcf_indels} ;
        fi
    """


rule bedtools_coverage_histogram:
  input:
    bam=lambda wildcards: f"data/{sample_to_reference[wildcards.sampleID]}/{wildcards.sampleID}/{wildcards.sampleID}.bam",
  output:
    histogram="4-bed_files/{sampleID}_genome_coverage_hist.tsv.gz",
  threads: 16
  conda:
    "envs/bedtools.yaml"
  shell:
    """
        mkdir -p 4-bed_files 4-bed_files/tmp
        tmpdir=$(mktemp -d 4-bed_files/tmp/{wildcards.sampleID}.hist.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        samtools view -b -F 4 -q 30 {input.bam} > "$tmpdir/filtered.bam"
        bedtools genomecov -ibam "$tmpdir/filtered.bam" | gzip -c > {output.histogram}
    """


rule bedtools_zero_coverage:
  input:
    bam=lambda wildcards: f"data/{sample_to_reference[wildcards.sampleID]}/{wildcards.sampleID}/{wildcards.sampleID}.bam",
  output:
    zero_coverage="4-bed_files/{sampleID}_merged_zero_covg_regions.tsv.gz",
  params:
    merge_distance=500,
  threads: 16
  conda:
    "envs/bedtools.yaml"
  shell:
    """
        mkdir -p 4-bed_files 4-bed_files/tmp
        tmpdir=$(mktemp -d 4-bed_files/tmp/{wildcards.sampleID}.zero.XXXXXX)
        trap 'rm -rf "$tmpdir"' EXIT
        samtools view -b -F 4 -q 30 {input.bam} > "$tmpdir/filtered.bam"
        bedtools genomecov -bga -ibam "$tmpdir/filtered.bam" > "$tmpdir/all_positions.tsv"
        awk -F'\t' '{{if ($NF == 0) print}}' "$tmpdir/all_positions.tsv" > "$tmpdir/zero_coverage.bed"
        if [ -s "$tmpdir/zero_coverage.bed" ]; then
            bedtools merge -d {params.merge_distance} -i "$tmpdir/zero_coverage.bed" > "$tmpdir/merged_zero_coverage.bed"
            bedtools coverage -b "$tmpdir/filtered.bam" -a "$tmpdir/merged_zero_coverage.bed" -hist |
                awk -F'\t' '{{if ($4 == 0) print}}' | gzip -c > {output.zero_coverage}
        else
            gzip -c /dev/null > {output.zero_coverage}
        fi
    """


rule mpileup2vcf_ancient:
  input:
    bam="data/{reference}/{sampleID}/{sampleID}.bam",
    ref="/nexus/posix0/MPIIB-keylab/reference_genomes/{reference}/genome.fasta",
  output:
    pileup="1-vcf/{sampleID}_ref_{reference}_aligned.sorted.pileup",
    variants="1-vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.vcf.gz",
    vcf_strain="1-vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.vcf.gz",
  group:
    'pileup_and_filter', 
  params:
    vcf_raw="1-vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.gz",
    minMAF = minMAF,
    regions = "{reference}_regions.bed",
  conda:
    "envs/samtools15_bcftools12.yaml",
  shell:
    """
        if [ -e {params.regions} ]; then 
            echo "regions file specified, generating mutations only on regions defined in regions.bed."
            samtools mpileup -l {params.regions} -q30 -x -s -O -d3000 -f {input.ref} {input.bam} > {output.pileup} ;
            samtools mpileup -l {params.regions} -q30 -t SP -d3000 -vf {input.ref} {input.bam} > {params.vcf_raw} ;
            bcftools call -c -Oz -o {output.vcf_strain} {params.vcf_raw} ;
            bcftools view -Oz -v snps -q {params.minMAF} {output.vcf_strain} > {output.variants} ;
            tabix -p vcf {output.variants} ;
            rm {params.vcf_raw}
        else 
            echo "regions.bed file does not exist, output will be across full genome";
            samtools mpileup -q30 -x -s -O -d3000 -f {input.ref} {input.bam} > {output.pileup} ;
            samtools mpileup -q30 -t SP -d3000 -vf {input.ref} {input.bam} > {params.vcf_raw} ;
            bcftools call -c -Oz -o {output.vcf_strain} {params.vcf_raw} ;
            bcftools view -Oz -v snps -q {params.minMAF} {output.vcf_strain} > {output.variants} ;
            tabix -p vcf {output.variants} ;
            rm {params.vcf_raw}
        fi
    """

# strain.vcf ==> vcf_to_quals.m ==>.quals.npz
rule vcf2quals_ancient:
  input:
    vcf_strain = rules.mpileup2vcf_ancient.output.vcf_strain,
  params:
    refGenomeDir="/nexus/posix0/MPIIB-keylab/reference_genomes/{reference}/"
  output:
    file_quals = "2-quals/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.quals.npz",
  group:
    'pileup_and_filter',
  run:
    from vcf_to_quals_snakemake import vcf_to_quals_snakemake
    vcf_to_quals_snakemake(sample_path_to_vcf = input.vcf_strain, sample_path_to_quals = output.file_quals, REF_GENOME_DIRECTORY = params.refGenomeDir)

# strain.pileup ==> pileup_to_diversity.m ==> diversity.mat
rule pileup2diversity_matrix_ancient:
  input:
    pileup = rules.mpileup2vcf_ancient.output.pileup,
  params:
    refGenomeDir="/nexus/posix0/MPIIB-keylab/reference_genomes/{reference}/",
  output:
    file_diversity = "3-diversity/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.diversity.npz",
    file_coverage = "3-diversity/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.coverage.npz",
  group:
    'pileup_and_filter',
  run:
    from pileup_to_diversity_matrix_snakemake import pileup_to_div_matrix_snakemake
    pileup_to_div_matrix_snakemake(sample_path_to_pileup = input.pileup, sample_path_to_diversity =  output.file_diversity, sample_path_to_coverage = output.file_coverage, ref_genome_directory = params.refGenomeDir)

rule remove_pileup_ancient:
  input:
    pileup = rules.mpileup2vcf_ancient.output.pileup,
  params:
    file_diversity = "3-diversity/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.diversity.npz",
    file_coverage = "3-diversity/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.coverage.npz",
  group:
    'pileup_and_filter',
  shell:
    "rm {input.pileup};"

rule cleanUp_ancient:
  input:
    part1 = expand("2-quals/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.quals.npz", zip, sampleID=SAMPLE_ls, reference=REF_Genome_ls),  # input not used, only required so snakemake waits with clean up until the end
    part2 = expand("3-diversity/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.diversity.npz", zip, sampleID=SAMPLE_ls, reference=REF_Genome_ls),  # input not used, only required so snakemake waits with clean up until the end
  output:
    "cleanUp_done.txt",
  shell:
    " touch {output} ;"

rule generate_next_samplescsv:
  input: 
    csv = "samples.csv",
  output:
    case_csv = "samples_case.csv",
  shell: 
    """ echo 'Path,Sample,ReferenceGenome,Outgroup' > {output.case_csv} ;"""
    " dir=$(pwd) ;"
    """ awk -v dir="$dir" 'BEGIN{{FS=OFS=","}} NR>1 {{print dir,$2,$3,$5}}' {input.csv} >> {output.case_csv} ;"""
