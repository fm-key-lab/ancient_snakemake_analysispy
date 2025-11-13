#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import gzip
from scipy.stats import ttest_ind, fisher_exact,ttest_1samp
from math import log10, floor
from analysispy_module import genomestats


## input output name parsing
def get_input_output_names(sample_path_to_pileup,sample_path_to_diversity,sample_path_to_coverage):
    fname_in=os.getcwd() +  '/' + sample_path_to_pileup  # Aro_Change: new; print this
    # Path to diversity.mat (output)
    fname_out=os.getcwd() +  '/' + sample_path_to_diversity  # Aro_Change: new; print this
    # Path to coverage.mat (secondary output)
    fname_out_cov=os.getcwd() +  '/' + sample_path_to_coverage  # Aro_Change: new; print this
    #fname_in='strain.pileup'; # Aro_Change: old
    #fname_out='diversity.mat'; # Aro_Change: old
    return (fname_in,fname_out,fname_out_cov)

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return floor(n*multiplier + 0.5) / multiplier    

## main function to be called in run block in snakemake:
def pileup_to_div_matrix_snakemake(sample_path_to_pileup,sample_path_to_diversity,sample_path_to_coverage,ref_genome_directory, min_reads_on_strand=20): ## default value in tami's code: 
    # Declare constants
    Phred_offset=33; #mpileup is told that the reads are in fastq format and corrects that so its always at Phred+33 when the mpileup comes out
    nts=[ord(nt) for nt in 'ATCGatcg']
    indels=[ord(i) for i in '+-']
    numfields=40; #CHANGED 2019.01
    numtests=50; #how many tests to do in parallel ##NOTE: CURRENTLY NOT IMPLEMENTED
    indelregion=3
    ## get input and output filenames to parse
    (fname_in,fname_out,fname_out_cov)=get_input_output_names(sample_path_to_pileup,sample_path_to_diversity,sample_path_to_coverage)
    # Get reference genome information
    #load alignment_info % Aro_Change: old
    #refgenome=ae.Genome; % Aro_Change: old    
    print(f'\n Reference genome: {ref_genome_directory} \n')
    [ChrStarts,GenomeLength,ScafNames] = genomestats(ref_genome_directory) # Aro_Change: new
    print(f'\n {GenomeLength} \n') # Aro_Change: print this
    #[ChrStarts,GenomeLength,ChromosomeIndicator,ScafNames]=genomestats(['/scratch/mit_lieberman/Reference_Genomes/'
    #refgenome]); %RUN_ON_CLUSTER % Aro_Change: old
    # %initialize
    ## maybe grab Genome length from give genomestats object?
    data=np.zeros((numfields,GenomeLength),dtype=int)
    line_dict={}
    with open(fname_in,'rt') as fid:
        index=0
        for line in fid:
            temp=np.zeros((numfields,1),dtype=int)
            if type(line) == str:
                parsed=line.strip().split('\t') ## split into list and get rid of any EOL 
                chrom=parsed[0]
                if len(ChrStarts) == 1:
                    position=int(parsed[1])-1 ## python is 0 based, not 1 based as in matlab
                else:
                    if chrom not in ScafNames:
                        raise Exception("scaffold name in pileup not found in reference")
                    position=int(ChrStarts[np.where(ScafNames == chrom)] + int(parsed[1])) - 1 ## first part finds start position of this chrom and adds loc of SNP on chrom // MF: changed from .index to np.where and changed int(line[1]) to int(parsed[1])
                if ord(parsed[2]) in nts: ## check first to make sure ref is ATCG or atcg, not N, Y, (other ambiguous calls)
                    ref=nts.index(ord(parsed[2])) ## get numerical value of reference base
                    if ref > 3:
                        ref = ref-4
                else: 
                    ref = None ##may want to be more permissive for parsing these, currently all base calls will not be counted for this position
                calls_ascii_values=np.array([ord(i) for i in parsed[4]])
                bq_ascii_values=np.array([ord(i) for i in parsed[5]]) #base quality, BAQ corrected
                mq_ascii_values=np.array([ord(i) for i in parsed[6]]) #mapping quality
                ## distance from tail, needs loop due to weird encoding of certain chars
                if parsed[7] == "*":
                    p_7=[0]
                else: 
                    p_7=parsed[7].split(',')
                td_ascii_values=np.array([int(i) for i in p_7]) #distance from tail 

                #starts of reads
                start_of_read_vector=np.where(calls_ascii_values==94)[0] ## 94 == '^' == start of read segment ## needs to find all positions
                for i in start_of_read_vector: #(do this to the ^ char and the next)
                    calls_ascii_values[i]=-1 ## convert start of read to -1
                    calls_ascii_values[i+1]=-1 #remove mapping character (mapping qual= ASCII in index after ^ - 33),
                    #absolutely required because the next chracter could be $
                
                #read ends ##TODO: CHECK THIS WITH A LINE WITH END READS
                end_of_read_vector=np.where(calls_ascii_values==36)[0] ## 36 == '$' == end of read segment 
                temp[37]=len(start_of_read_vector) + len(end_of_read_vector)
                calls_ascii_values[end_of_read_vector]=-1 ## convert end of reads to -1

                ## indels and calls from reads supporting indels 
                ## TODO: check indel vector where these appear
                indel_vector=np.concatenate([np.where(calls_ascii_values==43)[0], np.where(calls_ascii_values==45)[0]]) ## all positions that called an indel
                for i in indel_vector:
                    # find size 
                    forward_looking_index_for_indel_size = 1
                    while calls_ascii_values[i + forward_looking_index_for_indel_size] >=48 and calls_ascii_values[i + forward_looking_index_for_indel_size] <58:
                        forward_looking_index_for_indel_size+=1
                    if forward_looking_index_for_indel_size > 1:
                        indel_ascii=list(calls_ascii_values[i+1:i+forward_looking_index_for_indel_size])
                        indel_size=int(''.join(map(chr,indel_ascii)))
                        indeld=forward_looking_index_for_indel_size-1
                    else:
                        indel_size=calls_ascii_values[i+1]-48
                        indeld=1
                    # record that indel was found in nearby region
                    #store directly into data, as it affects lines earlier and later
                    if calls_ascii_values[i]==45: ## deletion
                        if (position - indelregion + 1) > 0 and (position+indel_size+indelregion) <= GenomeLength:
                            data[38:40,position-indelregion+1:position+indel_size+indelregion + 1 ]=data[38:40,position-indelregion+1:position+indel_size+indelregion+1]+1
                        elif (position - indelregion+1) > 0:
                            data[38:40,position-indelregion+1:]=data[38:40,position-indelregion+1:]+1
                        else:
                            data[38:40,0:position+indelregion+indel_size+1]=data[38:40,0:position+indel_size+indelregion+1]+1
                    else: ## insertion 
                        if (position-indelregion+1) > 0 and (position+indelregion) <= GenomeLength:
                            data[38,position-indelregion+1:position+indelregion+1]=data[38,position-indelregion+1:position+indelregion+1]+1
                        elif (position-indelregion+1)>0:
                            data[38,position-indelregion+1:]=data[38,position-indelregion+1:]+1
                        else:
                            data[38,0:position+indel_size+1]=data[38,0:position+indel_size+1]+1
                    #remove these calls from counting
                    calls_ascii_values[i:(i+indeld+indel_size+1)]=-1 # dont remove base that preceeds and indel
                #   ignore '*', as these deletion markers won't be counted towards score
                #   %qualities, etc and the presence of a deletion is recorded by the
                #   upstream '-' indicator 
                ## replace reference with actual calls
                if ref != None:
                    calls_ascii_values[np.where(calls_ascii_values == 46)[0]] = nts[ref] ## '.' forward, store in nts index values ATCGatcg==0:7
                    calls_ascii_values[np.where(calls_ascii_values == 44)[0]] = nts[ref+4] ## ',' reverse, store in nts index values ATCGatcg==0:7         
                ## indx reads for finding scores ## TODO might be bugged
                simple_calls=calls_ascii_values[np.where(calls_ascii_values>0)]## indices of all non-filtered calls (see above)
                for i in range(len(nts)):
                    current_nt=nts[i]
                    current_nt_indices=np.where(simple_calls==current_nt)[0]
                    size_for_this_nt=len(current_nt_indices) ##TODO CHECK IF RIGHT WITH OTHER FILES
                    if size_for_this_nt != 0:
                        temp[i]=size_for_this_nt
                        temp[i+8]=round_half_up(sum(bq_ascii_values[(current_nt_indices)])/size_for_this_nt) - Phred_offset ## average phred score for calls of given base
                        temp[i+16]=round_half_up(sum(mq_ascii_values[(current_nt_indices)])/size_for_this_nt) - 33 ## average mapping qual for this given base
                        temp[i+24]=round_half_up(sum(td_ascii_values[(current_nt_indices)])/size_for_this_nt) ## average tail distance
                ## The following section is critical for metagenomic samples 
                # find major and nextmajor allele
                allele_index_summed_calls=[(x,temp[x]+temp[x+4]) for x in range(4)] # add T and t calls, C and c calls, etc. keep index of major allele as first index in tuple
                allele_index_summed_calls.sort(key=lambda x: x[1][0]) ## sort by number of calls for a given base
                n1 = allele_index_summed_calls[-1][0] # major allele (most calls) (returns actual index of base in nts)
                n2 = allele_index_summed_calls[-2][0] # next most abundant allele (returns actual index of base in nts)
                
                x=np.logical_or(simple_calls==nts[n1], simple_calls==nts[n1+4]) # boolean array of reads with major alelle
                y=np.logical_or(simple_calls==nts[n2], simple_calls==nts[n2+4]) # boolean array of reads with minor alelle
                
                if sum(temp[0:4]) > min_reads_on_strand and sum(temp[4:8])> min_reads_on_strand and sum(y)*200>sum(x): # only calcualte p values of there are greater than 20 reads on each strand and MAF < .995
                    bttests_x=bq_ascii_values[x]-Phred_offset
                    bttests_y=bq_ascii_values[y]-Phred_offset # gets base qualities for major/minor alleles
                    mttests_x=mq_ascii_values[x]-Phred_offset
                    mttests_y=mq_ascii_values[y]-Phred_offset # gets mapping qualities for major/minor alleles
                    fttests_x=td_ascii_values[(np.where(simple_calls==nts[n1])[0])]
                    fttests_y=td_ascii_values[(np.where(simple_calls==nts[n2])[0])] # gets tail distances from fwd reads for major/minor alleles
                    rttests_x=td_ascii_values[(np.where(simple_calls==nts[n1+4])[0])]
                    rttests_y=td_ascii_values[(np.where(simple_calls==nts[n2+4])[0])] # gets tail distances from rev reads for major/minor alleles
                    
                    ## NOTE: matlab ttest will output result if len(y) == 1, treating it as the hypothesized mean of the normal distribution compared to x 
                    bp=ttest_ind(bttests_x,bttests_y)[1]
                    mp=ttest_ind(mttests_x,mttests_y)[1]
                    fp=ttest_ind(fttests_x,fttests_y)[1]
                    rp=ttest_ind(rttests_x,rttests_y)[1]
                    # report pval of above t-tests as -log10
                        ## coverting values of 0 and nan for output as -1:
                    if bp == 0 or np.isnan(bp):
                        temp[33]=-1
                    else: temp[33]=round_half_up(-log10(bp)) 

                    if mp == 0 or np.isnan(mp):
                        temp[34]=-1
                    else: temp[34]=round_half_up(-log10(mp))

                    if fp == 0 or np.isnan(fp):
                        temp[35]=-1
                    else: temp[35]=round_half_up(-log10(fp))
                    
                    if rp == 0 or np.isnan(rp):
                        temp[36]=-1
                    else: temp[36]=round_half_up(-log10(rp))
                    
                    ## currently not broken but testing against matlab script fails since matlab script *is* broken.
                    p=fisher_exact(np.array([[temp[n1][0], temp[n2][0]],[temp[n1+4][0],temp[n2+4][0]]]), alternative='two-sided')[1] # fisher exact test for strand bias (contingency table = # of reads that are fwd/rev and support major/minor allele)
                    if p == 0:
                        temp[32]=-1
                    else:
                        temp[32]=round_half_up(-log10(p))
                ## store the data
                data[0:38,position]=np.ravel(temp[0:38]) 
                ## NOTE: column 38 should possibly be overwritten here, if following matlab script
                ## however, this seems weird and like a bug in their code.
                ## if we want to use indels, we should figure this out thoroughly
            index+=1
    coverage=sum(data[0:8,:])
    np.savez_compressed(fname_out, data = data.astype(int))
    np.savez_compressed(fname_out_cov, coverage = coverage.astype(int))
