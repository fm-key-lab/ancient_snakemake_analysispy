# 
# =============================================================================
## General script for parsing raw CMT to slightly procesesd CMT
#  Necessary for large CMTs  
# author - Ian Light-Maka
# =============================================================================
# SNP-based analyses - QC File:
# =============================================================================
# Parse and filter raw output from Candidate Mutation Table generation(Snakemake)
# # =============================================================================

## import libraries not gotten from analysispy_module
import argparse 
import sys
import os
import re
import pickle
import numpy as np
import gzip
import json


sys.path.append("./local_analysis_initial_qc")
import local_analysis_initial_qc_modules as apy


# parser

parser = argparse.ArgumentParser(
                    prog='local_analysis_initial_qc_main',
                    description='Runs initial QC to reduce the size of candidate mutation table, based on various filtering metrics')

parser.add_argument('-p', '--parameter_json', 
                    help='Parameters for running initial QC')

# helper functions
def summary_sample_check(coverage_all,filter_parameter_sample_across_sites):
    # Within sample checks, across all positions
    goodsamples_output_path=f'goodsamples.npz'
    if os.path.exists(f'{goodsamples_output_path}'):
        print(f'summary_sample_check: existing output file for this check! loading {goodsamples_output_path}')
        goodsamples = np.load(f'{goodsamples_output_path}')['arr_0']
    else:
        failed_mean_coverage =  np.mean(coverage_all, axis=0) < filter_parameter_sample_across_sites['min_average_coverage_to_include_sample'] 
        failed_median_coverage =  np.median(coverage_all,axis=0) < filter_parameter_sample_across_sites['min_median_coverage_to_include_sample']

        goodsamples=failed_median_coverage | failed_mean_coverage
        np.savez_compressed(f'{goodsamples_output_path}',goodsamples)

    return goodsamples

#TODO: add median coverage filter
def within_sample_checks(quals,maf,coverage_forward_strand,coverage_reverse_strand,indels,coverage,filter_parameter_site_per_sample,optional_within_sample_checks=[]):
    # Witin sample checks, pos by pos
    ## Filter per mutation
    failed_within_sample_output_path=f'failed_within_sample.npz'
    if os.path.exists(f'{failed_within_sample_output_path}'):
        print(f'within_sample_checks: existing output file for this check! loading {failed_within_sample_output_path}')
        failed_within_sample = np.load(f'{failed_within_sample_output_path}')['arr_0']
    else:
        failed_quals = (quals < filter_parameter_site_per_sample['min_qual_for_call'])
        failed_maf=(maf < filter_parameter_site_per_sample['min_maf_for_call'])
        failed_forward=(coverage_forward_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
        failed_reverse=(coverage_reverse_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
        failed_cov=(coverage_reverse_strand + coverage_forward_strand < filter_parameter_site_per_sample['min_cov_on_pos'])
        failed_indels=(indels > (0.5*coverage) )
        # summarize and output
        failed_within_sample = (failed_quals | failed_maf | failed_forward | failed_reverse | failed_cov | failed_indels)
        np.savez_compressed(f'{failed_within_sample_output_path}',failed_within_sample)
    return failed_within_sample

def recombinant_check(optional_filtering,p,mutantAF,ingroup_bool,ancient_bool,filter_parameter_site_across_samples,failed_any_QC):
    failed_recombinant_output_path=f'failed_recombinant.npz'
    if os.path.exists(f'{failed_recombinant_output_path}'):
        print(f'recombinant_check: existing output file for this check! loading {failed_recombinant_output_path}')
        failed_recombinants = np.load(f'{failed_recombinant_output_path}')['arr_0']
    else:
        recombination_distance=filter_parameter_site_across_samples['distance_threshold_recombination']
        recombination_correlation=filter_parameter_site_across_samples['correlation_threshold_recombination']
        if optional_filtering['recombination'] == 'Ancient':
            failed_recombinants = apy.findrecombinantSNPs(p,mutantAF[ : , ancient_bool ],recombination_distance,recombination_correlation, failed_any_QC[ : , ancient_bool ] )[1]
        elif optional_filtering['recombination'] == 'All':
            failed_recombinants = apy.findrecombinantSNPs(p,mutantAF[ : , ingroup_bool ],recombination_distance,recombination_correlation, failed_any_QC[ : , ingroup_bool ] )[1]
        else:
            failed_recombinants=np.zeros(p.shape)
        np.savez_compressed(f'{failed_recombinant_output_path}',failed_recombinants)
    return failed_recombinants

def metagenomic_checks(optional_filtering,p,calls,minorAF,ancient_bool):
    """
    TODO: implement in snakemake
    # Metagenomic checks
    # aDNA checks (coverage percentile, and )
    # generated by file /Users/ad_loris/Nextcloud/keylab/projects/il_bronzeage_pestis_evo/scripts_master/raw_data_processing/generate_bedcoverage_files.sh
    bed_histogram_path='bed_files/final_out_files/*_genome_coverage_hist.tsv.gz'
    bed_zero_covg_path = 'bed_files/final_out_files/*_merged_zero_covg_regions.tsv.gz'

    failed_genomic_islands = apy.filter_bed_0_cov_regions(bed_zero_covg_path,p,scafNames,chrStarts,sampleNames,filter_parameter_site_per_sample['max_prop_0_covg_ancient'])
    failed_genomic_islands[:,~ancient_bool]=False

    failed_coverage_percentile=apy.filter_bed_cov_hist(bed_histogram_path,p,scafNames,chrStarts,sampleNames,coverage,filter_parameter_site_per_sample['max_percentile_cov_ancient'],two_tailed=False,upper=True)
    failed_coverage_percentile[:,~ancient_bool]=False
    """
    failed_metagenomic_output_path='failed_metagenomic.npz'
    if os.path.exists(f'{failed_metagenomic_output_path}'):
        print(f'metagenomic_check: existing output file for this check! loading {failed_metagenomic_output_path}')
        failed_metagenomic = np.load(failed_metagenomic_output_path)['arr_0']
    else:
        ## Removing SNP calls that are very close to nearby heterozygosity per sample (10bp default)
        if optional_filtering['heterozygosity'] == 'All':
            failed_heterozygosity=apy.find_calls_near_heterozygous_sites(p,minorAF,10,0.1)
        elif optional_filtering['heterozygosity'] == 'Ancient':
            failed_heterozygosity=apy.find_calls_near_heterozygous_sites(p,minorAF,10,0.1)
            failed_heterozygosity[:,~ancient_bool]=False
        else: 
            failed_heterozygosity = np.full(calls.shape, False)
        failed_metagenomic = failed_heterozygosity
        np.savez_compressed(f'{failed_metagenomic_output_path}',failed_metagenomic)
        #TODO: add union of failed covg percentile and failed genomic islands
        # failed_heterozygosity = (failed_genomic_islands | failed_coverage_percentile | failed_heterozygosity)
    return failed_metagenomic

def identify_indels(indel_depth,indel_support,indel_filtering_params,indel_index_for_identites,indel_identities):
    # apply filtering parameters on indels
    indel_total_depth=np.nansum(indel_depth,axis=2)

    ## getting goodpos for indels based on freebayes calls
    ## indel_depth == (pos x samples x (reads supporting highest GL indel x reads supporting all other calls))
    has_indel_call_support=(indel_support > indel_filtering_params['min_gl_diff']) \
        & (indel_total_depth >= indel_filtering_params['min_covg_on_pos']) \
        & (indel_depth[:,:,0]/indel_total_depth > indel_filtering_params['major_allele_freq_indel_call']) \
        
    has_indel_call_support[(np.sum(np.isnan(indel_index_for_identites),axis=1)/indel_index_for_identites.shape[1] > indel_filtering_params['max_fraction_ambiguous_samples']),:] = False

    has_indel = has_indel_call_support & ~np.isnan(indel_index_for_identites) & (indel_index_for_identites !=0)

    indel_sizes_called = np.empty(has_indel.shape)
    for index,row in enumerate(indel_index_for_identites):
        ref_len = len(indel_identities[index][0])
        for sample_index,value in enumerate(row):
            if not np.isnan(value):
                indel_sizes_called[index,sample_index]=len(indel_identities[index][int(value)]) - ref_len
            else: 
                indel_sizes_called[index,sample_index]=0


    ## define goodpos indels by finding all positions which have >1 indel type segregating (indel(s) or indel(s) and 0)
    goodpos_indels=np.where(np.sum(has_indel,axis=1)>0)[0]

    candidate_indels=[]
    for i in goodpos_indels:
        if len(np.unique(indel_sizes_called[i][~np.isnan(indel_sizes_called[i])])) > 1:
            candidate_indels.append(i)
    return goodpos_indels,candidate_indels,indel_sizes_called

def site_filter_check(calls,optional_filtering,failed_recombinants,minorAF):
    failed_site_filter_output_path='failed_site_filt.npz'
    if os.path.exists(f'{failed_site_filter_output_path}'):
        print(f'site_filter_check: existing file exists! loading {failed_site_filter_output_path}')
        failed_optional_siteFilt = np.load(failed_site_filter_output_path)['arr_0']
    else:
        failed_optional_siteFilt = np.full(calls.shape, False)
        if optional_filtering['recombination'] != 'None':
            failed_optional_siteFilt[failed_recombinants ,:] = 4
        if optional_filtering['minoraf_covariance']:
            cov_scores=[]
            for x in minorAF:
                cov_scores.append(np.cov(x))
            minorAF_covariance_failed=np.where(np.array(cov_scores) > np.percentile(np.array(cov_scores),optional_filtering['minoraf_covariance']))[0]
            failed_optional_siteFilt[minorAF_covariance_failed,:]=4
        np.savez_compressed(f'{failed_site_filter_output_path}',failed_optional_siteFilt)
    return failed_optional_siteFilt

def blast_masking():
    pass

def generate_fasta(goodpos_final2useTree,calls,sampleNames,refnt,refgenome,analysis_params_output_name,name_append=''):
    # get data and filter for goodpos_final
    calls_for_treei = calls[ goodpos_final2useTree, : ]; 

    # build sampleNames  (that passed filter, see above from sampleNames_all --> sampleNames) w/ metainfo
    treesampleNamesLong = sampleNames # include all samples 

    # translate index to nucleotide
    calls_for_tree = apy.idx2nts(calls_for_treei) # ATCGN translation

    # add reference nucleotide for all positions
    refgenome_nts_for_tree=refnt[goodpos_final2useTree]
    calls_for_tree_ref_outgroup = np.concatenate((apy.idx2nts(refgenome_nts_for_tree[:, None]),calls_for_tree),axis=1)

    treesampleNamesLong_ref_outgroup = np.append([f'{refgenome}'],treesampleNamesLong)
    apy.write_calls_sampleName_to_fasta(calls_for_tree_ref_outgroup,treesampleNamesLong_ref_outgroup,f'{analysis_params_output_name}{name_append}')

def save_qc_filtered(goodpos_final,counts,quals,coverage_forward_strand,coverage_reverse_strand,refnti_m,p,refgenome,sampleNames,outgroup_bool,contig_positions,mutantAF,maf,maNT,minorNT,minorAF,calls,hasmutation,analysis_params_output_name):
    # output fully reduced and filtered CMT
    cmtFile_sub_gp = f'final_cmts/{analysis_params_output_name}_post_initialqc_candidate_mutation_table.pickle.gz'
    os.makedirs('final_cmts', exist_ok=True)
    with gzip.open(cmtFile_sub_gp, 'wb') as pickle_file:
        pickle.dump({
            'counts': counts[:, :, goodpos_final],
            'quals': quals[goodpos_final, :],
            'coverage_forward_strand': coverage_forward_strand[goodpos_final, :],
            'coverage_reverse_strand': coverage_reverse_strand[goodpos_final, :],
            'refnti_m': refnti_m[goodpos_final, :],
            'p': p[goodpos_final],
            'refgenome': refgenome,  
            'sampleNames': sampleNames,  
            'outgroup_bool': outgroup_bool,  
            'contig_positions': contig_positions[goodpos_final, :],
            'mutantAF': mutantAF[goodpos_final, :],
            'maf': maf[goodpos_final, :],
            'maNT': maNT[goodpos_final, :],
            'minorNT': minorNT[goodpos_final, :],
            'minorAF': minorAF[goodpos_final, :],
            'calls': calls[goodpos_final, :],
            'hasmutation': hasmutation[goodpos_final, :],
            'goodpos_raw_cmt': goodpos_final,
            'goodpos_cleaned_cmt': np.array(range(len(goodpos_final))),
            'analysis_params_output_name': analysis_params_output_name
            }, pickle_file)

    # quickly test that dumping worked correctly:
    with gzip.open(cmtFile_sub_gp, 'rb') as pickle_file:
        pickle_dictionary = pickle.load(pickle_file)
        counts_loaded = pickle_dictionary['counts']
        quals_loaded = pickle_dictionary['quals' ]
        coverage_forward_strand_loaded = pickle_dictionary['coverage_forward_strand' ]
        coverage_reverse_strand_loaded = pickle_dictionary['coverage_reverse_strand' ]
        refnti_m_loaded = pickle_dictionary['refnti_m' ]
        p_loaded = pickle_dictionary['p' ]
        refgenome_loaded = pickle_dictionary['refgenome' ]
        sampleNames_loaded = pickle_dictionary['sampleNames' ]
        outgroup_bool_loaded = pickle_dictionary['outgroup_bool' ]
        contig_positions_loaded = pickle_dictionary['contig_positions' ]
        mutantAF_loaded = pickle_dictionary['mutantAF' ]
        maf_loaded = pickle_dictionary['maf' ]
        maNT_loaded = pickle_dictionary['maNT' ]
        minorNT_loaded = pickle_dictionary['minorNT' ]
        minorAF_loaded = pickle_dictionary['minorAF' ]
        calls_loaded = pickle_dictionary['calls' ]
        hasmutation_loaded = pickle_dictionary['hasmutation' ]
        goodpos_final_raw_cmt_loaded = pickle_dictionary['goodpos_raw_cmt' ]
        goodpos_final_cleaned_cmt_loaded = pickle_dictionary['goodpos_cleaned_cmt' ]
        analysis_params_output_name_loaded = pickle_dictionary['analysis_params_output_name']


    print('Checking parsing for counts - - - - - - - - - - - - -  PARSING OK?:',np.all(counts[:, :, goodpos_final] == counts_loaded))
    print('Checking parsing for quals - - - - - - - - - - - - - - PARSING OK?:',np.all(quals[goodpos_final, :] == quals_loaded))
    print('Checking parsing for coverage_forward_strand - - - - - PARSING OK?:',np.all(coverage_forward_strand[goodpos_final, :] == coverage_forward_strand_loaded))
    print('Checking parsing for coverage_reverse_strand - - - - - PARSING OK?:',np.all(coverage_reverse_strand[goodpos_final, :] == coverage_reverse_strand_loaded))
    print('Checking parsing for refnti_m - - - - - - - - - - - -  PARSING OK?:',np.all(refnti_m[goodpos_final, :] == refnti_m_loaded))
    print('Checking parsing for p - - - - - - - - - - - - - - - - PARSING OK?:',np.all(p[goodpos_final] == p_loaded))
    print('Checking parsing for refgenome - - - - - - - - - - - - PARSING OK?:',np.all(refgenome  == refgenome_loaded))
    print('Checking parsing for sampleNames - - - - - - - - - - - PARSING OK?:',np.all(sampleNames  == sampleNames_loaded))
    print('Checking parsing for outgroup_bool - - - - - - - - - - PARSING OK?:',np.all(outgroup_bool  == outgroup_bool_loaded))
    print('Checking parsing for contig_positions - - - - - - - -  PARSING OK?:',np.all(contig_positions[goodpos_final, :] == contig_positions_loaded))
    print('Checking parsing for mutantAF - - - - - - - - - - - -  PARSING OK?:',np.all(mutantAF[goodpos_final, :] == mutantAF_loaded))
    print('Checking parsing for maf - - - - - - - - - - - - - - - PARSING OK?:',np.all(maf[goodpos_final, :] == maf_loaded))
    print('Checking parsing for maNT - - - - - - - - - - - - - -  PARSING OK?:',np.all(maNT[goodpos_final, :] == maNT_loaded))
    print('Checking parsing for minorNT - - - - - - - - - - - - - PARSING OK?:',np.all(minorNT[goodpos_final, :] == minorNT_loaded))
    print('Checking parsing for minorAF - - - - - - - - - - - - - PARSING OK?:',np.all(minorAF[goodpos_final, :] == minorAF_loaded))
    print('Checking parsing for calls - - - - - - - - - - - - - - PARSING OK?:',np.all(calls[goodpos_final, :] == calls_loaded))
    print('Checking parsing for hasmutation - - - - - - - - - - - PARSING OK?:',np.all(hasmutation[goodpos_final, :] == hasmutation_loaded))
    print('Checking parsing for goodpos_final_raw_cmt - - - - - - PARSING OK?:',np.all(goodpos_final == goodpos_final_raw_cmt_loaded))
    print('Checking parsing for goodpos_final_cleaned_cmt - - - - PARSING OK?:',np.all(np.array(range(len(goodpos_final))) == goodpos_final_cleaned_cmt_loaded))
    print('Checking parsing for analysis_params_output_name - - - PARSING OK?:', analysis_params_output_name == analysis_params_output_name_loaded)

def validate_json(pared_json):
    pass

def main(parameter_json):
    with open(f'{parameter_json}') as f_in:
        json_parsed=json.load(f_in)

    # validate all paths, etc, error if anything missing and return useful info to end user
    validate_json(json_parsed)

    # setup I/O
    analysis_params_output_name=json_parsed['input_output']['runtime_name']

    input_output_dir=json_parsed['input_output']['input_directory']
    
    os.makedirs(f'{input_output_dir}/{analysis_params_output_name}',exist_ok=True)
    os.chdir(f'{input_output_dir}/{analysis_params_output_name}')

    ref_genome_dir=json_parsed['input_output']['ref_directory']

    # Filtering parameters
    # =============================================================================
    filtering_dicts=json_parsed['filtering']

    # Much of this will need to vary according to your reference genome,
    # coverage, and particular samples
    # adjust within json!

    filter_parameter_sample_across_sites = filtering_dicts['filter_parameter_sample_across_sites']

    filter_parameter_site_per_sample = filtering_dicts['filter_parameter_site_per_sample']

    filter_parameter_site_across_samples = filtering_dicts['filter_parameter_site_across_samples']

    indel_filtering_params = filtering_dicts['indel_filtering_params']

    optional_filtering = filtering_dicts['optional_filtering']

    ######################################################
    ### SETUP DONE ###### SETUP DONE ###### SETUP DONE ###
    ######################################################
    # run full QC script

    ######################################################
    ## START PROCESSING DATA #### START PROCESSING DATA ##
    ######################################################
    # 
    # Read in genome information
    # =============================================================================
    [chrStarts, genomeLength, scafNames] = apy.genomestats(ref_genome_dir)
    refgenome=ref_genome_dir.split('/')[-1]

    # Load data from candidate_mutation_
    # =============================================================================
    [quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats,indel_p,indel_depth,indel_support,indel_identities,indel_index_for_identites] = apy.read_candidate_mutation_table_pickle_gzip(f'{input_output_dir}/candidate_mutation_table.pickle.gz')

    # 
    # Mean cov per samples based on all positions in counts (aka p)
    # =============================================================================

    mean_cov_p_per_sample = np.vstack( (np.array(['sample','mean_p_cov']),np.column_stack( (sampleNames,np.mean( np.sum(counts,axis=1),axis=1)) )) )    
    np.savetxt('./mean_cov_per_sample.csv', mean_cov_p_per_sample, delimiter=',',fmt='%s')

    median_cov_p_per_sample = np.vstack( (np.array(['sample','median_p_cov']),np.column_stack( (sampleNames,np.median( np.sum(counts,axis=1),axis=1)) )) )    
    np.savetxt('./median_cov_per_sample.csv', median_cov_p_per_sample, delimiter=',',fmt='%s')

    # =============================================================================
    # convert outputs to more permanent datatypes
    # =============================================================================

    sampleNames_all=np.asarray(sampleNames,dtype=object)
    quals_all=-quals;
    counts_all=counts;
    coverage_all = counts_all.sum(axis=1).transpose() # axis=1 == rows; transpose needed > rows: pos and col: samples
    indel_depth_all=indel_depth
    indel_support_all=indel_support

    # convert indels from pileups
    # =============================================================================
    # % indel counter reduction to indel count
    # indel_counter: The first statistic is the number of reads (at this position and in this sample) that 
    # support an indel. The second statistics is the number of reads (at this position and in this sample) that support a deletion.
    # for now we need only the count for indels overall 
    indel_counter = indel_counter[:,0,:].transpose()
    # indel counter >> 50% indel row:1; row2 >> deletion

    # Display samples the NOT fullfill min_average_coverage_to_include_sample:
    # =============================================================================
    #
    # ID samples of outgroup/ingroup 
    if len(json_parsed['sample_labelling']['outgroup_name']) > 0:
        outgroup_sample = json_parsed['sample_labelling']['outgroup_name']
        outgroup_bool = np.isin(sampleNames,[f'{outgroup_sample}'])
    elif len(json_parsed['sample_labelling']['outgroup_pattern']) > 0:
        outgroup_pattern = re.compile(fr'{json_parsed["sample_labelling"]["outgroup_pattern"].replace}')
        outgroup_bool = np.isin(sampleNames,list(filter(outgroup_pattern.match, sampleNames)))
    else:
        outgroup_bool=np.zeros(len(sampleNames)).astype(bool)

    ancient_pattern = re.compile(r'^SP\.*|^M219')
    ancient_bool = np.isin(sampleNames,list(filter(ancient_pattern.match, sampleNames)))


    # =============================================================================
    # Define goodsamples and filter data, good samples have > avg coverage
    # =============================================================================

    goodsamples = summary_sample_check(coverage_all,filter_parameter_sample_across_sites)

    #Breakpoint: Too few samples passed filter, checking that at least 2 samples pass QC
    if np.sum(goodsamples) < 2:
        print("Too few samples fullfill filter criteria! >> skip: " + refgenome)

    sampleNames = sampleNames_all[goodsamples]
    counts = counts_all[goodsamples , : , : ] # keep only level (samples) that fullfil filter!
    quals = quals_all[ : , goodsamples ]
    coverage = coverage_all[ : ,goodsamples]
    indels = indel_counter[:,goodsamples]
    outgroup_bool = outgroup_bool[goodsamples]
    ancient_bool = ancient_bool[goodsamples]


    ## TODO: update how this is done
    # TODO: CAN ALSO JUST GRAB OUTGROUP BY OUTPUTTING THE SAMPLES CSV from case
    outgroup_idx=np.nonzero(outgroup_bool)[0]

    # ingroup array (bool, idx) used later
    ingroup_bool = np.invert(outgroup_bool)
    ingroup_idx = np.nonzero(ingroup_bool)[0]


    indel_depth=indel_depth_all[:,goodsamples,:]
    indel_support=indel_support_all[:,goodsamples]
    indel_index_for_identites=indel_index_for_identites[:,goodsamples]
    indel_total_depth=np.nansum(indel_depth,axis=2)

    num_samples = len(sampleNames)

    coverage_forward_strand = counts[:,0:4,:].sum(axis=1).transpose()
    coverage_reverse_strand = counts[:,4:8,:].sum(axis=1).transpose()

    indel_support[:,outgroup_idx]=np.nan
    #indel_sizes_called[:,outgroup_idx]=np.nan

    #goodpos_indels,candidate_indels,indel_sizes_called=identify_indels(indel_depth,indel_support,indel_filtering_params,indel_index_for_identites,indel_identities)
    
    # =============================================================================
    # Extract refnt and define out/in-group bools
    # =============================================================================

    ## Note ancnti/outs_nti defined below after filtered calls has been generated!
    ## get reference allele for all p; NOTE: analysis.m stored ATCG as double-digit-numeric
    # use ref nt for mutation calls. important if multiple outgroups called 
    refnt = apy.extract_outgroup_mutation_positions(ref_genome_dir, apy.p2chrpos(p,chrStarts));
    refnti = apy.nts2idx(refnt)
    refnti_m = np.tile(refnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele



    # Find positions with fixed mutations
    # 
    # Extract allele and frequencies
    contig_positions = apy.p2chrpos(p,chrStarts) # 1col: chr, 2col: pos on chr; for all p
    [maf, maNT, minorNT, minorAF] = apy.div_major_allele_freq(counts) 
    calls = maNT
    intermediate_hasmutation=(calls != refnti_m) & (calls < 4)
    intermediate_hasmutation[:,outgroup_bool] = False
    print("Initial number of candidate variants: ", len(np.where( np.sum(intermediate_hasmutation, axis=1) > 0 )[0]))

    # NOTE: function assumes first 8 rows in counts == 4nucl fwd&rev! watch out if extended counts used!
    # NOTE: maf==0 -> no data;minorAF==0 -> no minor allele/ or no major allele; NT number corresponds to index in NTs [ATCG] or if maf==0 > NA == 4  
        
    #  Make some basic structures for finding mutations
    mutantAF = np.zeros(maNT.shape)
    mutantAF[maNT != refnti_m] = maf[ maNT != refnti_m]; 


    # Define mutations we do not trust in each and across samples.
    # goodpos are indices of p that we trust

    failed_within_sample=within_sample_checks(quals,maf,coverage_forward_strand,coverage_reverse_strand,indels,coverage,filter_parameter_site_per_sample)

    calls[failed_within_sample]=4

    intermediate_hasmutation_post_sample_checks=(calls != refnti_m) & (calls < 4)
    intermediate_hasmutation_post_sample_checks[:,outgroup_bool] = False
    intermediate_num_variants=len(np.where( np.sum(intermediate_hasmutation_post_sample_checks, axis=1) > 0 )[0])
    variants_removed_failed_within_sample=len(p)-intermediate_num_variants
    print("Number of variants after this filtering step: ", intermediate_num_variants)
    print("Number of variants removed: ", variants_removed_failed_within_sample)

    failed_metagenomic=metagenomic_checks(optional_filtering,p,calls,minorAF,ancient_bool)
    calls[(failed_metagenomic | failed_within_sample)]=4

    intermediate_hasmutation_post_sample_checks_and_metagenomics=(calls != refnti_m) & (calls < 4)
    intermediate_hasmutation_post_sample_checks_and_metagenomics[:,outgroup_bool] = False
    intermediate_num_variants=len(np.where( np.sum(intermediate_hasmutation_post_sample_checks_and_metagenomics, axis=1) > 0 )[0])
    variants_removed_failed_within_sample_failed_metagenomics=len(p)-variants_removed_failed_within_sample-intermediate_num_variants

    print("Number of variants after this filtering step: ", intermediate_num_variants)
    print("Number of variants removed: ", variants_removed_failed_within_sample_failed_metagenomics)

    # save positions and samples where a metagenomic variant was filtered
    np.where(intermediate_hasmutation_post_sample_checks_and_metagenomics != intermediate_hasmutation_post_sample_checks)

    with open(f'{analysis_params_output_name}_masked_metagenomic_snvs.tsv', 'w') as f:
        f.write(f'Masked metagenomics snvs on ancient samples')
        f.write('\n')
        for x,y in zip(apy.p2chrpos(p[np.where(intermediate_hasmutation_post_sample_checks_and_metagenomics != intermediate_hasmutation_post_sample_checks)[0]],chrStarts),sampleNames[np.where(intermediate_hasmutation_post_sample_checks_and_metagenomics != intermediate_hasmutation_post_sample_checks)[1]]):
            true_chrom=scafNames[x[0]-1]
            position_on_chrom=x[1]+1
            f.write(f'{y}\t{true_chrom}:{position_on_chrom}\n')
    

    # Across sample checks 
    ## Remove putative recombinants
    failed_any_QC = ( failed_metagenomic | failed_within_sample)

    failed_recombinants=recombinant_check(optional_filtering,p,mutantAF,ingroup_bool,ancient_bool,filter_parameter_site_across_samples,failed_any_QC)


    ## Filter per site across samples
    # Ignore here outgroup samples!
    siteFilt = np.any(( (calls[:,ingroup_bool]>3).sum(axis=1) >= ((num_samples-np.sum(outgroup_bool)) * filter_parameter_site_across_samples['max_fraction_ambiguous_samples']) \
                            ,np.median( coverage[:,ingroup_bool], axis=1) < filter_parameter_site_across_samples['min_median_coverage_position'] ),axis=0)
    calls[ siteFilt ,:] = 4 # sites that fail qc -> 4, for all samples incl. outgroup     


    failed_optional_siteFilt = site_filter_check(calls,optional_filtering,failed_recombinants,minorAF)
    
    calls[ failed_optional_siteFilt ] = 4 # sites that fail qc -> 4, for all samples incl. outgroup     


    # NOTE: func below takes forever with many SNPs...saved below
    # TODO: update mutation quality to check if path exists
    [mutQual, mutQualIsolates] = apy.ana_mutation_quality(calls[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
    mutQual = np.nan_to_num(mutQual, nan=-1) # turn mutQual nan's to -1; necessary to avoid later warning

    # translate filtered calls of ingroup into goodpos. mutations we believe. fixedmutation part removed in v6.
    hasmutation = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!
    hasmutation[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
    candpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

    # calculate last set of removed variants
    variants_removed_across_sites=intermediate_num_variants-len(candpos)
    print("Number of variants after this filtering step: ", len(candpos))
    print("Number of variants removed: ", variants_removed_across_sites)

    goodpos = candpos
    print(goodpos.size,'goodpos found.')

    if json_parsed['input_output']['save_fasta']:
        generate_fasta(goodpos,calls,sampleNames,refnt,refgenome,analysis_params_output_name,name_append='')


    if optional_filtering['run_blast_masking']:
        blast_masking()


    save_qc_filtered(goodpos,counts,quals,coverage_forward_strand,coverage_reverse_strand,refnti_m,p,refgenome,sampleNames,outgroup_bool,contig_positions,mutantAF,maf,maNT,minorNT,minorAF,calls,hasmutation,analysis_params_output_name)

    """
    #
    # # # # 
    # # # # # # # # # # # # # # # # 
    # ancient singleton validation # 
    # # # # # # # # # # # # # # # # 
    # # # #
    ancient_singletons = np.intersect1d(candpos,np.where((np.sum(hasmutation[:,np.in1d(sampleNames,ancient_sample_names)],axis=1)==1) & (np.sum(hasmutation[:,~np.in1d(sampleNames,ancient_sample_names)],axis=1)==0))[0])

    for pos in ancient_singletons:
        samples=[x for x in sampleNames[hasmutation[p==p[pos]][0]] if x in ancient_sample_names]
        print(p[pos], samples)
        

    sample_pos_to_output={}
    for pos in ancient_singletons:
        samples=[x for x in sampleNames[hasmutation[p==p[pos]][0]] if x in ancient_sample_names][0]
        if samples in sample_pos_to_output:
            sample_pos_to_output[samples].append(pos)
        else:
            sample_pos_to_output[samples] = [pos]

    regions_to_inspect=[]
    for s in sample_pos_to_output:
        for singletons_position in sample_pos_to_output[s]:
            chrom,position_on_chrom=apy.p2chrpos([p[singletons_position]],chrStarts)[0]
            true_chrom=scafNames[chrom-1]
            true_position_on_chrom_1_based=position_on_chrom+1
            regions_to_inspect.append(f'{s},{true_chrom}:{true_position_on_chrom_1_based}-{true_position_on_chrom_1_based}')
    np.savetxt('singletons_to_extract_mappinginfo.txt', np.array(regions_to_inspect,dtype=str), fmt='%s')


    # BREAKPOINT FOR RUNNING BLAST+
    # TODO: integrate
    # run blast on identified ancient singleton reads
    # local_analysis/analysis_py/blast_execution_ancient_singletons.sh

    # # # # # # # # # # # # # # # # # # # # # # # # #

    # 

"""


## execute 

args = parser.parse_args()
main(args.parameter_json)

"""
# LOAD BACK IN RESULTS
def reparse_singleton_mpileup_calls_to_reads(blast_result_file):
    # TODO: move to apy
    fields=['Chrom','Pos','Ref','Covg','Calls','BQ','MQ','Read_Pos','Reads']
    pileup_file=pd.read_csv(blast_result_file.replace('_out.tsv.gz','.pileup').replace('blast_results','mpileup_with_readnames'),sep='\t', names=fields)
    reads=pileup_file.Reads.iloc[0].split(',')
    reference_call=pileup_file.Ref.iloc[0]
    # reparse calls
    raw_calls=pileup_file.Calls.iloc[0].replace('.',reference_call).replace(',',reference_call.lower())
    cleaned_calls=[]
    indices_to_skip=set()
    for idx,call in enumerate(raw_calls):
        if call in ['+', '-']: # skip indexing for all indel info
            size_of_indel=raw_calls[idx+1]
            indices_to_skip.add([i for i in range(idx,idx+size_of_indel+2)]) # exclude indices for eg +,2,A,A
        if call in ['^','$']: # skip indexing for all start/end info
            indices_to_skip.add(idx)
        if idx not in indices_to_skip:
            cleaned_calls.append(call)
    if len(cleaned_calls) != len(reads):
        raise ValueError
    return np.array(cleaned_calls),np.array(reads)

def update_counts_for_singleton(counts,call_support_remaining,sample_index_this_query,p_index_this_query):
    # update counts, calls, maf, and cov with cleaned read support for ancient singletons
    nt_order=['A','T','C','G','a','t','c','g']
    for idx,nt in enumerate(nt_order):
        cleaned_counts_this_nt=np.sum(call_support_remaining==nt)
        counts[sample_index_this_query,idx,p_index_this_query]=cleaned_counts_this_nt
    return counts


def update_matrices_post_singleton_counts_update(counts_updated,calls_previous,maNT_previous,refnti_m,filter_parameter_site_per_sample,sample_pos_to_output):
    # Recalculate entire datastructures    
    coverage_forward_strand_updated = counts_updated[:,0:4,:].sum(axis=1).transpose()
    coverage_reverse_strand_updated = counts_updated[:,4:8,:].sum(axis=1).transpose()
    [maf_updated, maNT_updated, minorNT_updated, minorAF_updated] = apy.div_major_allele_freq(counts_updated) 
    calls_updated = maNT_updated
    mutantAF_updated = np.zeros(maNT_updated.shape)
    mutantAF_updated[maNT_updated != refnti_m] = maf_updated[ maNT_updated != refnti_m]; 
    failed_maf_updated=(maf_updated < filter_parameter_site_per_sample['min_maf_for_call'])
    failed_forward_updated=(coverage_forward_strand_updated < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
    failed_reverse_updated=(coverage_reverse_strand_updated < filter_parameter_site_per_sample['min_cov_per_strand_for_call'])
    failed_cov_updated=(coverage_reverse_strand_updated + coverage_forward_strand_updated < filter_parameter_site_per_sample['min_cov_on_pos'])
    failed_any_QC_updated = ( failed_maf_updated | failed_forward_updated | failed_reverse_updated | failed_cov_updated )
    calls_updated[failed_any_QC_updated] = 4
    # get indices to subset to only ancient singleton positions investigated:
    sample_indices_to_update_from_dict=[]
    pos_indices_to_update_from_dict=[]
    for sample in sample_pos_to_output:
        sample_index=np.where(sampleNames==sample)[0][0]
        sample_indices_to_update_from_dict+=[sample_index]*len(sample_pos_to_output[sample])
        pos_indices_to_update_from_dict+=sample_pos_to_output[sample]
    sample_indices_to_update_from_dict=np.array(sample_indices_to_update_from_dict)
    pos_indices_to_update_from_dict=np.array(pos_indices_to_update_from_dict)

    # only update to calls, maNT for the positions that are investigated
    calls_previous [ pos_indices_to_update_from_dict,sample_indices_to_update_from_dict ] = calls_updated[pos_indices_to_update_from_dict,sample_indices_to_update_from_dict]
    
    return [calls_previous,coverage_forward_strand_updated,coverage_reverse_strand_updated,maf_updated, calls_previous, minorNT_updated, minorAF_updated, mutantAF_updated]



    blast_results=glob.glob('ancient_singleton_validation/blast_results/*.tsv.gz')

    pseudo_taxids=[1649845,633,109458,273123,349747,502800,502801,748676,748679,1218101,1286084,1286085,1286086,1286087,1286088,1286089,1324587,1324608,1355443,748672]
    pestis_taxids=[1035377,1234662,1345701,1345702,1345703,1345704,1345705,1345706,1345707,1345708,1345710,1455696,187410,214092,229193,349746,360102,377628,386656,547048,632,637382,637386,649716,748678]

    singletons_to_mask=[]
    singletons_to_mask_possibly_fine=[]
    singletons_to_mask_sample=[]
    counts_should_update=0
    reads_failed={}
    for blast_result_file in blast_results:
        blast_results_dict={}
        sampleid_for_query,chromid_posid_for_query=blast_result_file.split('/')[2].split('_NC_')
        chromid='NC_'+chromid_posid_for_query.split(':')[0]
        query=int(chromid_posid_for_query.split(':')[1].split('-')[0])
        cleaned_calls,cleaned_reads=reparse_singleton_mpileup_calls_to_reads(blast_result_file)
        pileup_file=pd.read_csv(blast_result_file.replace('_out.tsv.gz','.pileup').replace('blast_results','mpileup_with_readnames'),sep='\t', names=['Chrom','Pos','Ref','Covg','Calls','BQ','MQ','Read_Pos','Reads'])
        pileup_coverage=pileup_file.Covg.iloc[0]
        reads_in_counts=0
        with gzip.open(blast_result_file, 'rt') as f:
            for l in f:
                query_id,hit_id,hit_taxid,bitscore,pident=l.strip().split('\t')
                if query_id in pileup_file.Reads.iloc[0].split(','):
                    reads_in_counts+=1
                    if query_id not in blast_results_dict:
                        blast_results_dict[query_id]={(bitscore,pident):set([int(hit_taxid)])}
                    elif (bitscore,pident) in blast_results_dict[query_id]:
                        blast_results_dict[query_id][(bitscore,pident)].add(int(hit_taxid))
        # parse what the MAF is in the supported read call
        p_index_this_query=np.where(p==apy.chrpos2p(chromid,int(query)-1,scafNames,chrStarts))[0][0]
        if p_index_this_query not in ancient_singletons:
            print('indexing error')
        sample_index_this_query=np.where(sampleNames==sampleid_for_query)[0][0]
        maf_this_call=maf[p_index_this_query,sample_index_this_query]
        # number of pestis/pseudotb specific reads MUST exceed this proportion of reads to be retained
        reads_passed=[]
        reads_failed[sampleid_for_query]=[]
        reads_assessed=0
        for read_assessed in blast_results_dict:
            reads_assessed+=1
            for q,taxid_hits in blast_results_dict[read_assessed].items():
                taxid_hits_array=np.array([x for x in taxid_hits])
                if len(np.intersect1d(np.array(pseudo_taxids+pestis_taxids),taxid_hits_array)) > 0:
                    reads_passed.append(read_assessed)
                else:
                    reads_failed[sampleid_for_query].append(read_assessed)
        num_reads_passed=len(reads_passed)
        counts_coverage=np.sum(counts[sample_index_this_query,:,p_index_this_query])
        if num_reads_passed != counts_coverage:
            counts_should_update+=1
        print(f'{sampleid_for_query} {chromid}:{query} -- Passed: {num_reads_passed}, Assessed: {reads_assessed}, Counts Coverage: {counts_coverage}')
        if 0 < num_reads_passed < pileup_coverage:
            singletons_to_mask_possibly_fine.append(p_index_this_query)
            #print(f'Passed: {num_reads_passed}, Assessed: {reads_assessed}, Total Coverage: {pileup_coverage}')
        call_support_remaining=cleaned_calls[np.in1d(cleaned_reads,np.array(reads_passed))]
        counts=update_counts_for_singleton(counts,call_support_remaining,sample_index_this_query,p_index_this_query)

    # outputting reads that failed check of origin to pestis or pseudotb
    with open('blast_failed_reads_with_sample.tsv','w') as f:
        for s in reads_failed:
            if len(reads_failed[s])>0:
                reads_failed_parsed="\t".join(reads_failed[s])
                f.write(f'{s}\t{reads_failed_parsed}\n')

    with open('blast_failed_reads.tsv','w') as f:
        for s in reads_failed:
            if len(reads_failed[s])>0:
                for r in reads_failed[s]:
                    f.write(f'{r}\n')

    # finding sample + position indices that have been updated
    samples_indices_to_update=np.where(hasmutation[np.unique(np.where(counts!=counts_copy)[2])])[1]
    indices_to_update=np.unique(np.where(counts!=counts_copy)[2])
    # NOTE: ONLY these positions should be possibly updated in following function call

    #
    # Now final update to relevant matrices for saving
    [calls, 
    coverage_forward_strand, 
    coverage_reverse_strand, 
    maf, 
    maNT,
    minorNT, 
    minorAF, 
    mutantAF ] = update_matrices_post_singleton_counts_update(counts,calls,maNT,refnti_m,filter_parameter_site_per_sample,sample_pos_to_output)

    # confirm no positions other than those set above have changed:
    print("Only changed positions are within the set of indices to update?:",len(np.where(calls[indices_to_update,samples_indices_to_update]!=calls_copy[indices_to_update,samples_indices_to_update])[0]) == len(np.where(calls!=calls_copy)[0]))

    hasmutation_final = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!

    hasmutation_final[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
    candpos_final = np.where( np.sum(hasmutation_final, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

    goodpos_final = candpos_final
    print('Number of ancient singletons removed: ',len(goodpos)-len(goodpos_final))
    # saving indices of discarded ancient singletons:
    with open(f'{analysis_params_output_name}_masked_singletons.tsv', 'w') as f:
        f.write(f'Masked singletons on ancient samples')
        f.write('\n')
        for x,y in zip(apy.p2chrpos(p[np.where(calls!=calls_copy)[0]],chrStarts),sampleNames[np.where(calls!=calls_copy)[1]]):
            true_chrom=scafNames[x[0]-1]
            position_on_chrom=x[1]+1
            f.write(f'{y}\t{true_chrom}:{position_on_chrom}\n')


    print(goodpos_final.size,'goodpos found.')

    print("Carry out interactive QC investigations (below) for these positions and adjust parameters as necssary.")
"""
