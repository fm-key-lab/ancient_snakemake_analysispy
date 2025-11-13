# -*- coding: utf-8 -*-
"""
---Gathers everything together for candidate_mutation_table---
NOTE: Still reads in many *.mat files etc. Further purging of matlab necessary!

Output:
# path_candidate_mutation_table: where to write
# candidate_mutation_table.mat, ex. results/candidate_mutation_table.mat

---


# Inputs (changed to argparse usage):
     path_to_p_file: where to find all_positions.mat
     path_to_sample_names_file: where to find text file with sample names
         (space delimited)
     path_to_outgroup_boolean_file: where to find text file with outgroup
         booleans (space delimited, 1=outgroup, 0=not)
    path_to_list_of_quals_files: where to find text file with list of
       quals.mat files for each sample (space delimited)
     path_to_list_of_diversity_files: where to find text file with list of
       diversity.mat files for each sample (space delimited)
# Output:
     path_candidate_mutation_table: where to write
     candidate_mutation_table.mat, ex. results/candidate_mutation_table.mat

# Note: All paths should be relative to pwd!


## Version history

     This is adapted from TDL's build_mutation_table_master_smaller_file_size_backup.m
  #   Arolyn, 2018.12.19: This script was written as part of the transition to snakemake. 
          It performs the part of the case step that gathers data for
         Quals and counts and saves candidate_mutation_table.mat
  #   Arolyn, 2019.02.12: Added another matlab variable that stores indel
          statistics called 'indel_counter'.
  #   Tami, 2019.12.12: Converted into python and also added ability save coverage data
  #   Felix: 2020.01-04: Continous Debugged and adapted script for streamlined Snakemake implementation. 
  #                      Added argparse for proper argument parsing and optional coverage matrix build.
"""

''' load libraries '''
from itertools import count
import numpy as np
import pickle
import scipy.io as sio
import os
import sys,argparse
import gzip
from scipy import sparse
import vcf
from Bio import SeqIO
import glob

''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='''\
                            Gathers everything together for candidate_mutation_table.
                            Optional: Builds coverage matrix (optional w/ double-standardized matrix)
                               ''',
                               epilog="Questions or comments? --> fkey@mit.edu")
parser.add_argument("-r", dest="refgenomedir", help="Reference genome folder",required=True,action='store')
parser.add_argument("-p", dest="allpositions", help="All positions p file (*mat)",required=True,action='store')
parser.add_argument("-s", dest="sampleNames", help="File with sample names",required=True,action='store')
parser.add_argument("-g", dest="outgroupBool", help="String outgroup bool",required=True,action='store')
parser.add_argument("-q", dest="qualfiles", help="String qual matrix paths",required=True,action='store')
parser.add_argument("-d", dest="divfiles", help="String diversity paths",required=True,action='store')
parser.add_argument("-i", dest="indelvcffiles", help="String indel vcf path",required=True,action='store')
parser.add_argument("-o", dest="candidate_mutation_table", help="Output candidate mutation table. Py pickle structure (*.pickle.gz)",required=True,action='store')
parser.add_argument("-c", dest="get_cov", help="Set flag to build raw coverage matrix as sparse csr gzip numpy object (dirname+cov_raw_sparsecsr_mat.npz)",action="store_true", default=False)
parser.add_argument("-n", dest="get_dbl_norm_cov", help="Set flag to build double normalized coverage matrix as sparse csr gzip numpy object (dirname+cov_norm_sparsecsr_mat.npz)",action="store_true", default=False)
args = parser.parse_args()


'''Functions'''

def genomestats(REFGENOMEFOLDER):
    # parse ref genome to extract relevant stats
    # accepts genome.fasta or genome.fasta.gz (gzip) in refgenomefolder
    fasta_file = glob.glob(REFGENOMEFOLDER + '/genome.fasta')
    if len(fasta_file) != 1:
        fasta_file_gz = glob.glob(REFGENOMEFOLDER + '/genome.fasta.gz')
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOMEFOLDER)
        else: # genome.fasta.gz
            refgenome = SeqIO.parse(gzip.open(fasta_file_gz[0], "rt"),'fasta')
    else: # genome.fasta
        refgenome = SeqIO.parse(fasta_file[0],'fasta')
    Genomelength = 0
    ChrStarts = []
    ScafNames = []
    for record in refgenome:
        ChrStarts.append(Genomelength) # chr1 starts at 0 in analysis.m
        Genomelength = Genomelength + len(record)
        ScafNames.append(record.id)
    # turn to np.arrys!
    ChrStarts = np.asarray(ChrStarts,dtype=int)
    Genomelength = np.asarray(Genomelength,dtype=int)
    ScafNames = np.asarray(ScafNames,dtype=object)
    return [ChrStarts,Genomelength,ScafNames]

def main(path_to_refgenome_dir,path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files, path_to_list_of_diversity_files, path_to_indel_vcf, path_to_candidate_mutation_table, flag_cov_raw_sparse_matrix,flag_cov_norm_sparse_scale_matrix):
    # get refgenome info:
    [ChrStarts,Genomelength,ScafNames]=genomestats(path_to_refgenome_dir)

    pwd=os.getcwd()
    
    # p: positions on genome that are candidate SNPs
    print('Processing candidate SNP positions...')
    
    
    infile=np.load(path_to_p_file) # from previous step, should include variable called p
    p=infile['p'].flatten()
    p=p-1 #since converting from MATLAB!!!
    print('Total number of positions: ' + str(len(p)))
    
    
    # SampleNames: list of names of all samples
    print('Processing sample names...')
    
    fname =  pwd + '/' + path_to_sample_names_file  
    fid = open( fname, "r" ) # Input is space separated text file, in one line
    SampleNames = fid.readline().split()
    SampleNames = np.array(SampleNames)
    fid.close()
    
    numSamples = len(SampleNames) # save number of samples
    print('Total number of samples: ' + str(numSamples))
    
    
    ## in_outgroup: booleans for whether or not each sample is in the outgroup
    print('Processing outgroup booleans...')
    
    fname =  pwd + '/' + path_to_outgroup_boolean_file  
    fid = open( fname, "r" ) # Input is space separated text file, in one line
    in_outgroup_string = fid.readline().split()
    in_outgroup=np.array(in_outgroup_string)
    in_outgroup = in_outgroup.reshape(1,len(in_outgroup)) # reshape 2d array for analysis.py: 1row and numSamples cols
    fid.close()
    outgroup_idx=np.where(in_outgroup == '1')[1]
    
    
    ## Quals: quality score (relating to sample purity) at each position for all samples
    print('Gathering quality scores at each candidate position...')
    # Import list of directories for where to quals.mat for each sample
    fname =  pwd + '/' + path_to_list_of_quals_files  
    fid = open( fname, "r" ) # Input is space separated text file, in one line
    paths_to_quals_files = fid.readline().split()
    fid.close()
     
    
    # Make Quals
    Quals = np.zeros((len(p), numSamples), dtype='int') # initialize
    for i in range (numSamples):
        print('Loading quals matrix for sample: ' + str(i)) 
        print('Filename: ' + paths_to_quals_files[i]) 
        infile=np.load(paths_to_quals_files[i]) # from previous step, should include variable called p
        quals=infile['quals'].flatten()
        Quals[:,i]=quals[p]
    
    
    
    ## counts: counts for each base from forward and reverse reads at each candidate position for all samples
    print('Gathering counts data at each candidate position...\n')
    
    # Import list of directories for where to diversity.mat for each sample
    fname =  pwd + '/' + path_to_list_of_diversity_files  
    fid = open( fname, "r" ) # Input is space separated text file, in one line
    paths_to_diversity_files = fid.readline().split()
    fid.close()
    
    tempfile=np.load(paths_to_diversity_files[1]) 
    data=tempfile['data']
    size=np.shape(data)
    GenomeLength=size[1]
        
    # Make counts and coverage at the same time
    counts = np.zeros((8, len(p), numSamples),dtype='uint') # initialize
    all_coverage_per_bp = np.zeros((numSamples,GenomeLength),dtype='uint') # Added 2019.12.12
    indel_counter=np.zeros((2, len(p), numSamples), dtype='uint') # Added 2019.02.12
    for i in range (numSamples):
        print('Loading counts matrix for sample: ' + str(i)) 
        print('Filename: '+ paths_to_diversity_files[i]) 
        infile=np.load(paths_to_diversity_files[i]) 
        data=infile['data']
        counts[:,:,i]=data[0:8,p]
        if flag_cov_raw_sparse_matrix:
            all_coverage_per_bp[i,:]=sum(data[0:8,:])
        indel_counter[:,:,i]=data[38:40,p] # Added 2019.02.12 reads supporting indels and reads supporting deletions
    counts = counts.transpose(2,0,1) # counts reshape for analysis.py: 0:samples,1:ACTG,2:p
    indel_counter = indel_counter.transpose(2,0,1) # indel_counter reshape for analysis.py: 0:samples,2,p 
    
    #print('Getting all the coverage information...\n')
    #[all_coverage_per_bp, ~, all_maf_per_bp] = get_all_coverage(SampleInfo, GenomeLength)
    
    # Normalize coverage by sample and then position; ignore /0 ; turn resulting inf to 0
    
    if flag_cov_norm_sparse_scale_matrix:
        with np.errstate(divide='ignore',invalid='ignore'):
            array_cov_norm = ( all_coverage_per_bp - np.mean(all_coverage_per_bp,axis=1,keepdims=True) ) / np.std(all_coverage_per_bp,axis=1,keepdims=True) # ,keepdims=True maintains 2D array (second dim == 1), necessary for braodcasting
            array_cov_norm[ ~np.isfinite(array_cov_norm) ] = 0
            
            # 2nd normalisation
            array_cov_norm = ( array_cov_norm - np.mean(array_cov_norm,axis=0,keepdims=True) ) / np.std(array_cov_norm,axis=0,keepdims=True) # ,keepdims=True maintains 2D array (second dim == 1), necessary for braodcasting
            array_cov_norm[ ~np.isfinite(array_cov_norm) ] = 0
    

    ## turn into sparse csr matrices for more efficient computation
    # scale norm matrix by 1000 and save as int64 to slim matrix as much as possible
    # save matrices
    if os.path.dirname(path_to_candidate_mutation_table) == '': # make sure cov matrix goes to same folder as cmt
        outdir = ''
    else:
        outdir = os.path.dirname(path_to_candidate_mutation_table) + '/'
    if flag_cov_raw_sparse_matrix:
        all_coverage_per_bp_csr = sparse.csr_matrix(all_coverage_per_bp)
        sparse.save_npz(outdir+'cov_raw_sparsecsr_mat.npz', all_coverage_per_bp_csr,compressed=True)
    if flag_cov_norm_sparse_scale_matrix:
        array_cov_norm_scaled_csr = sparse.csr_matrix((np.round(array_cov_norm,3)*1000).astype('int64'))
        sparse.save_npz(outdir+'cov_norm_sparsecsr_mat.npz', array_cov_norm_scaled_csr,compressed=True)
    
    ## NOTE: coverage_stats for all pos in diversity file (eg only positions with some coverage) has indexes 0-10 with num sites that have covg in that bin
    ## NOTE: coverage_stats for all pos in diversity file has index 11 with median coverage across all sites
    ## NOTE: coverage_stats for all pos in diversity file has index 12 with mean coverage across all sites
    ## NOTE: coverage stats for all pos in diversity file has index 13 with standard deviation of coverage across all sites
    ## NOTE: Due to this calculation, the coverage estimates are slightly elevated for mean and standard deviation
    coverage_stats=np.zeros((numSamples,14),dtype='uint') 
    for index in range(numSamples):
        unique,counts_for_dict=np.unique(all_coverage_per_bp[index], return_counts=True)
        counts_dict=dict(zip(unique,counts_for_dict))
        ## find mean and median, create bins 0-10 four outputting
        counts_for_output=[0]*11
        total_covg=0
        total_pos=0
        for val in counts_dict:
            if val <= 10:
                counts_for_output[val]=counts_dict[val]
            else:
                counts_for_output[10]+=counts_dict[val]
            total_pos+=counts_dict[val]
            total_covg+=val*counts_dict[val]
        mean_covg=total_covg/total_pos
        median_val=total_pos/2
        count_to_median=0
        index_for_median=0
        while count_to_median < median_val:
            count_to_median+=counts_for_dict[index_for_median]
            index_for_median+=1
        median=index_for_median-1
        ## save into coverage_stats
        coverage_stats[index,0:11]=counts_for_output
        coverage_stats[index,11]=median
        coverage_stats[index,12]=mean_covg
        coverage_stats[index,13]=np.std(all_coverage_per_bp[index])

    print("Parsing indels into CMT")
    print(f'Loading file',path_to_indel_vcf)
    ## Parse indels into CMT
    vcf_reader = vcf.Reader(filename=path_to_indel_vcf)
     ## make 2d matrices for outputting, p x samples
        # indel_p (p) - start position of variant called in joint indel vcf
        # depth (3d, pos x samples x (support alt, support non-alt) ) - support highest likelihood alt, support non-alt
        # indel_support (quals) - difference in genotype likelihood between (indel)variant and reference 
        # indel_size (calls) - size of variant (#bases added if +, #bases lost if -) for highest likelihood (indel)variant

    collector_indel_support=[]
    collector_indel_depth=[]
    collector_index_for_identity=[]
    collector_indel_identity=[]
    collector_indel_pos = []
    for record in vcf_reader:
        position_on_chr=record.POS
        chrom=record.CHROM
        position = ChrStarts[chrom == ScafNames] + position_on_chr - 1 # consistent 0indexing for CMT p matrix
        ref_len=len(record.REF)
        variant_impact=np.array([len(alt)-ref_len for alt in record.ALT])
        indices_indel=np.where(variant_impact != 0)[0]
        indel_identites=[record.REF]+[str(alt) for alt in record.ALT]
        if len(indices_indel) != 0: ## complex variants may not include an ins/del
            depth_at_pos=np.zeros((len(SampleNames),2)) 
            indel_support_at_pos=np.zeros((len(SampleNames)))
            indel_index_for_identity=np.zeros((len(SampleNames)))
            for sample_record in record.samples:
                sample_name=sample_record.sample
                sample_name_index=np.where(SampleNames==sample_name)
                if sample_record.data.GL: ## only parse sample data if there are reads to cover position
                    #get genotype likelihoods
                    genotype_likelihoods=np.array(sample_record.data.GL) 
                    #calculate difference in ML between top two likelihoods (closest to 0)
                    order_of_support=np.argsort(genotype_likelihoods)
                    ml_difference=genotype_likelihoods[order_of_support[-1]]-genotype_likelihoods[order_of_support[-2]]
                    depth_at_indel_all=sample_record.data.DP
                    reads_for_max_support=sample_record.data.AD[order_of_support[-1]] ## index +1, in AD (reads supporting each call) ref is first entry
                    ## fill corresponding CMT dataframes (will be converted to numpy arrays later)
                    depth_at_pos[sample_name_index]=np.array([reads_for_max_support,depth_at_indel_all-reads_for_max_support])
                    indel_support_at_pos[sample_name_index]=ml_difference
                    indel_index_for_identity[sample_name_index]=order_of_support[-1] # record what ref,alt1,alt2... index is highest support
                else:
                    depth_at_pos[sample_name_index]=np.nan
                    indel_support_at_pos[sample_name_index]=np.nan
                    indel_index_for_identity[sample_name_index]=np.nan
            # set outgroup call ( as np.nan)
            depth_at_pos[outgroup_idx]=np.nan
            indel_support_at_pos[outgroup_idx]=np.nan
            indel_index_for_identity[outgroup_idx]=np.nan
            # fill possible things, but only if there is an acutal supported alt call
            if not np.all(indel_index_for_identity[~np.isnan(indel_index_for_identity)]==0):
                collector_indel_support.append(indel_support_at_pos)
                collector_indel_depth.append(depth_at_pos)
                collector_index_for_identity.append(indel_index_for_identity)
                collector_indel_identity.append(indel_identites)
                collector_indel_pos.append(position)
    ## convert to numpy arrays for output
    
    indel_support=np.array(collector_indel_support,dtype='float')
    indel_depth=np.array(collector_indel_depth,dtype='float')
    indel_index_for_identity=np.array(collector_index_for_identity,dtype='float')
    indel_p=np.array(collector_indel_pos,dtype='int')
    indel_identites=np.array(collector_indel_identity,dtype='object')
    ## Save cmt!   
    with gzip.open(path_to_candidate_mutation_table, 'wb') as f: 
        pickle.dump([SampleNames, p, counts, Quals, in_outgroup, indel_counter, coverage_stats,indel_p,indel_depth,indel_support,indel_identites,indel_index_for_identity], f,protocol=4) # protocol=4 for storage of files >4gb
    
    print('DONE')


if __name__ == "__main__":
    path_to_refgenome_dir=args.refgenomedir
    path_to_p_file=args.allpositions
    path_to_sample_names_file=args.sampleNames
    path_to_outgroup_boolean_file=args.outgroupBool
    path_to_list_of_quals_files=args.qualfiles
    path_to_list_of_diversity_files=args.divfiles
    path_to_indel_vcf=args.indelvcffiles
    path_to_candidate_mutation_table=args.candidate_mutation_table
    flag_cov_raw_sparse_matrix=args.get_cov
    flag_cov_norm_sparse_scale_matrix=args.get_dbl_norm_cov
    if flag_cov_norm_sparse_scale_matrix and not flag_cov_raw_sparse_matrix:
        flag_cov_raw_sparse_matrix = True
        print('Selected to build double normalized coverage matrix. Raw coverage matrix will be build, too.')
    main(path_to_refgenome_dir,path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files, path_to_list_of_diversity_files, path_to_indel_vcf, path_to_candidate_mutation_table, flag_cov_raw_sparse_matrix,flag_cov_norm_sparse_scale_matrix)

