#########################################
# LIEBERMAN LAB SNAKEFILE FOR CASE STEP #
#########################################

# Version History:
# # 2022.03 Martin: Snakemake can now deal with multiple reference genomes at once
# # 2020.01 Felix: Added build_candidate_mutation_table.py incl. output of whole-genome coverage matrix and double-normalized matrix.
# # 2019.02.12, Arolyn: rule candidate_mutation_table now calling make_candidate_mutation_table_snakemake_indels.m instead of make_candidate_mutation_table_snakemake.m
# # 2018.12, Felix/Arolyn: Original Snakefile from the lab hackathon

# Reminders:
# # Put your own email address in cluster.slurm.json so that you get emails about failures. No one else wants your slurm emails.

import sys,os,subprocess
SCRIPTS_DIRECTORY = "./scripts"
sys.path.insert(0, SCRIPTS_DIRECTORY)

## USER defined variables (in theory do not need to be touched)
spls = "samples_case.csv"
maxFQ = -30 # threshold (from mapping quality) for including position

## pre-snakemake
from case_modules import *

''' PRE-SNAKEMAKE '''
## define couple of lists from samples_case.csv
## FORMAT: Path,Sample,ReferenceGenome,Outgroup
# NOTE: case analysis expects:
#   - unique sample IDs
#   - Col1==Paths, points to the snakemake directory used for mapping, not raw data! > only single entries accepted!
[PATH_ls, SAMPLE_ls, REF_Genome_ls, OUTGROUP_ls] = read_samples_caseCSV(spls)

## Te following lines are essential to run the snakemake with multiple reference genomes and get a candidate mutation table for each reference genome
## generate two dictories, which have as key the unique reference geneome for all samples and their outgroup tags 
SAMPLE_per_ref_dict = unique_key_generator(REF_Genome_ls, SAMPLE_ls)
OUTGROUP_per_ref_dict = unique_key_generator(REF_Genome_ls, OUTGROUP_ls)
REF_Genomes = list(set(REF_Genome_ls))

''' SNAKEMAKE '''
rule all:
    input:
        # # Data links only # #
        # expand("data/vcf/{sampleID}_ref_{reference}_outgroup{outgroup}.vcf.gz",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
        # expand("data/qual/{sampleID}_ref_{reference}_outgroup{outgroup}.quals.mat",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
        # expand("data/diversity/{sampleID}_ref_{reference}_outgroup{outgroup}.diversity.mat",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
        # # Everything # #
        expand("4-candidate_mutation_table/{reference}/candidate_mutation_table.pickle.gz",reference=set(REF_Genome_ls)) ## for each unique entry (set()) the output is expected
        # # Including cleanup step # #
        #"logs/DONE_cleanUp",

rule variants2positions:
    input:
        variants = "1-vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.vcf.gz",
    params:
        REF_GENOME_DIRECTORY = "/nexus/posix0/MPIIB-keylab/reference_genomes/{reference}/",
        outgroup_tag = "{outgroup}", # boolean (0==ingroup or 1==outgroup)
    output:
        mat_positions = "0-temp_pos/{sampleID}_ref_{reference}_outgroup{outgroup}_positions.npz",
    group:
        "var2pos",
    run:
        generate_positions_single_sample_snakemake_withoutbool(input.variants, output.mat_positions, maxFQ, params.REF_GENOME_DIRECTORY, params.outgroup_tag)


## combination steps are done per reference genome, double expand is therefore needed!
rule combine_positions_prep:
    input:
        mat_positions = lambda wildcards: expand(expand("0-temp_pos/{{sampleID}}_ref_{reference}_outgroup{{outgroup}}_positions.npz",reference=wildcards.reference),zip,sampleID=SAMPLE_per_ref_dict[wildcards.reference],outgroup=OUTGROUP_per_ref_dict[wildcards.reference]), ## Expand at first over the wildcard reference genome and mask the other two wildcarfs (done with the '{{}}'. Then, as just the zipped output is expected, zip the other two wildcards together
    params:
        outgroup_bool = lambda wildcards:  expand( "{outgroup}" , outgroup=OUTGROUP_per_ref_dict[wildcards.reference] ), ## expand for each reference genome individually and store in separate output files
    output:
        string_input_p_positions = "0-temp_pos/{reference}/string_file_other_p_to_consider.txt",
        string_outgroup_bool = "0-temp_pos/{reference}/string_outgroup_bool.txt",
    group:
        "pre_cand_mut_table",
    run:
        with open(output.string_input_p_positions, "w") as f:
            print(*input.mat_positions, sep=" ", file=f)
        with open(output.string_outgroup_bool, "w") as f:
            print(*params.outgroup_bool, sep=" ", file=f)


rule combine_positions:
    input:
        string_input_p_positions = "0-temp_pos/{reference}/string_file_other_p_to_consider.txt",
        string_outgroup_bool = "0-temp_pos/{reference}/string_outgroup_bool.txt",
    params:
        file_other_p_to_consider = "add_positions/other_positions.npz", ## candidate positions which could be added manually 
        REF_GENOME_DIRECTORY = "/nexus/posix0/MPIIB-keylab/reference_genomes/{reference}/"
    output:
        mat_positions = "0-temp_pos/{reference}/allpositions.npz",
    group:
        "pre_cand_mut_table",
    run:
        combine_positions_snakemake_out(input.string_input_p_positions, params.file_other_p_to_consider, output.mat_positions, input.string_outgroup_bool, params.REF_GENOME_DIRECTORY, maxFQ)


# build input for candidate_mutation_table
rule string_diversity_mat:
    input:
        cluster_in = rules.combine_positions.output.mat_positions, ## Needed to streamline rules and be able to group them on hpc cluster per ref genome
        diversity_mat = lambda wildcards: expand(expand("3-diversity/{{sampleID}}_ref_{reference}_aligned.sorted.strain.variant.diversity.npz",reference=wildcards.reference),sampleID=SAMPLE_per_ref_dict[wildcards.reference]),
        #diversity_mat = expand("data/diversity/{sampleID}_ref_{reference}_outgroup{outgroup}.diversity.npz",zip,sampleID=SAMPLE_ls, reference=REF_Genome_ls, outgroup=OUTGROUP_ls),
    output:
        string_diversity_mat = "0-temp_pos/{reference}/string_diversity_mat.txt",
    group:
        "pre_cand_mut_table",
    run:
        with open(output.string_diversity_mat, "w") as f:
            print(*input.diversity_mat, sep=" ", file=f)


# build input for candidate_mutation_table
rule string_quals_mat:
    input:
        cluster_in = rules.string_diversity_mat.output.string_diversity_mat, ## Needed to streamline rules and be able to group them on hpc cluster per ref genome
        quals_mat = lambda wildcards: expand(expand("2-quals/{{sampleID}}_ref_{reference}_aligned.sorted.strain.variant.quals.npz",reference=wildcards.reference),sampleID=SAMPLE_per_ref_dict[wildcards.reference]),
    output:
        string_qual_mat = "0-temp_pos/{reference}/string_qual_mat.txt",
    group:
        "pre_cand_mut_table",
    run:
        with open(output.string_qual_mat, "w") as f:
            print(*input.quals_mat, sep=" ", file=f)



# build input for candidate_mutation_table
## cannot be grouped as this is a single independent step
rule string_sampleID_names:
    input:
        rules.string_quals_mat.output.string_qual_mat, ## Needed to streamline rules and be able to group them on hpc cluster per ref genome
    params:
        sampleID_names = lambda wildcards: expand( "{sampleID}" , sampleID=SAMPLE_per_ref_dict[wildcards.reference] ), ## expand for each reference genome individually and store in separate output files
    output:
        string_sampleID_names = "0-temp_pos/{reference}/string_sampleID_names.txt",
    group:
        "pre_cand_mut_table",
    run:
        with open(output.string_sampleID_names, "w") as f:
            print(*params.sampleID_names, sep=" ", file=f)


rule candidate_mutation_table:
    input:
        mat_positions = "0-temp_pos/{reference}/allpositions.npz",
        string_diversity_mat = "0-temp_pos/{reference}/string_diversity_mat.txt",
        string_qual_mat = "0-temp_pos/{reference}/string_qual_mat.txt",
        string_sampleID_names = "0-temp_pos/{reference}/string_sampleID_names.txt",
        string_outgroup_bool = "0-temp_pos/{reference}/string_outgroup_bool.txt",
        string_indel_vcf = "1-vcf/ref_{reference}_non_outgroup_indels_complex.vcf.gz"
    output:
        candidate_mutation_table = "4-candidate_mutation_table/{reference}/candidate_mutation_table.pickle.gz",
        raw_sparse_cov_matrix = "4-candidate_mutation_table/{reference}/cov_raw_sparsecsr_mat.npz",
        double_norm_sparse_cov_matrix = "4-candidate_mutation_table/{reference}/cov_norm_sparsecsr_mat.npz", ## sparse coverage matrix normalized over samples and positions 
    params:
        REF_GENOME_DIRECTORY = "/nexus/posix0/MPIIB-keylab/reference_genomes/{reference}/"
    group:
        "cand_mut_table_clean",
    conda:
        "envs/py_for_cmt.yaml",
    shell:
        # -c/-n optional flag to build cov/norm matrix in folder of cmt. check -h for help.
        """
        python3 scripts/build_candidate_mutation_table_npz_script.py -r {params.REF_GENOME_DIRECTORY} -p {input.mat_positions} -s {input.string_sampleID_names} -g {input.string_outgroup_bool} -q {input.string_qual_mat} -d {input.string_diversity_mat} -i {input.string_indel_vcf} -o {output.candidate_mutation_table} -cn
        """


rule cleanUp:
    input:
        candidate_mutation_table = rules.candidate_mutation_table.output.candidate_mutation_table
    params:
        temp_folder = "0-temp_pos"
    output:
        "cleanUp_done.txt"
    group:
        "cand_mut_table_clean",
    shell:
        " rm -rf {params.temp_folder} ; touch {output} "
