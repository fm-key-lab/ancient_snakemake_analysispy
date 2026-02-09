
import gzip
import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import csv
import glob
import pandas as pd
import pickle
from collections import OrderedDict
from Bio import Phylo
import subprocess
from scipy import stats
import math
from matplotlib import rc
import matplotlib.pyplot as plt
from Bio.Data import CodonTable
import collections
from functools import partial # for gff_parse() and plot_coverage_fwd_rev_stacked()
from statsmodels.stats.proportion import proportion_confint
from statsmodels.stats.multitest import multipletests
import networkx


###################################################################################################
## # # # # # # # # # # # # # # # # # Import/Output of data # # # # # # # # # # # # # # # # # # # ##
###################################################################################################
# Import of Snakemake CMT
def read_candidate_mutation_table_pickle_gzip(file_cmt_pickle_gz):
    # read the candidate_mutation_table.pickle.gz files
    with gzip.open(file_cmt_pickle_gz, 'rb') as f:
        cmt = pickle.load(f)
        sampleNames = np.array(cmt[0])
        p = np.array(cmt[1]) # p 0-based index of genomic position (means actual position is p+1)
        counts = np.array(cmt[2])
        quals = np.array(cmt[3])
        in_outgroup = np.array(cmt[4])
        indel_counter = np.array(cmt[5])
        coverage_stats = np.array(cmt[6]) ## each row is sample, col [0-10) == covg bins 1x, 2x... >10x; [10]==median covg; [11]==mean; [12]==stddev 
        indel_p = np.array(cmt[7])
        indel_depth = np.array(cmt[8])
        indel_support = np.array(cmt[9])
        indel_identites = np.array(cmt[10])
        indel_index_for_identities = np.array(cmt[11])
    return [quals,p,counts,in_outgroup,sampleNames,indel_counter,coverage_stats,indel_p,indel_depth,indel_support,indel_identites,indel_index_for_identities]

# Import of genome info
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

# Data output (fasta generation)
def write_calls_sampleName_to_fasta(calls_for_tree,treeSampleNames,timestamp):
    fa_file = open(timestamp+".fa", "w")
    for i,name in enumerate(treeSampleNames):
        nucl_string = "".join(list(calls_for_tree[:,i]))
        fa_file.write(">" + name + "\n" + nucl_string + "\n")
    fa_file.close()


###################################################################################################
## # # # # # # # # # # Data manipulations (moving to and from datatypes) # # # # # # # # # # # # ##
###################################################################################################
def nts2idx(nts_array,missingdata="?"):
    # translate nucl array to array containing nucleotide indices
    nucl = np.array(['A','T','C','G',missingdata],dtype=object) # add 5th element --> no data! == index 4
    indices = [0,1,2,3,4] # values present in index-array
    for i in indices:
        nts_array[ nts_array == nucl[i] ] = i
    return nts_array.astype(int)

def idx2nts(calls,missingdata="?"):
    # translate index array to array containing nucleotides
    nucl = np.array(['A','T','C','G',missingdata],dtype=object) # add 5th element --> no data! == index 4
    palette = [0,1,2,3,4] # values present in index-array
    index = np.digitize(calls.ravel(), palette, right=True)
    return nucl[index].reshape(calls.shape)

def chrpos2p(chrom,pos,chrNames, ChrStarts):
    '''# return 2col array with chr and pos on chr
    #p...continous, ignores chr
    #pos: like p, 0-based'''
        
    # get chr and pos-on-chr
    this_chrom=np.where(chrNames==chrom)[0]
    p_to_return=ChrStarts[this_chrom]+pos
    return p_to_return

def p2chrpos(p, ChrStarts):
    '''# return 2col array with chr and pos on chr
    #p...continous, ignores chr
    #pos: like p, 0-based'''
        
    # get chr and pos-on-chr
    chr = np.ones(len(p),dtype=int)
    if len(ChrStarts) > 1:
        for i in ChrStarts[1:]:
            chr = chr + (p > i) # when (p > i) evaluates 'true' lead to plus 1 in summation. > bcs ChrStarts start with 0...genomestats()
        positions = p - ChrStarts[chr-1] # [chr-1] -1 due to 0based index
        pos = np.column_stack((chr,positions))
    else:
        pos = np.column_stack((chr,p))
    return pos

###################################################################################################
## # # # # # # # # # # # # # # # # Converting raw matrices # # # # # # # # # # # # # # # # # # # ##
###################################################################################################
def extract_outgroup_mutation_positions(REFGENOMEFOLDER,position,CMTpy=True):
    # extracts the ref nucleotide for every position. positions needs to be sorted by chr
    # reads genome.fasta or genome.fasta.gz    
    fasta_file = glob.glob(REFGENOMEFOLDER + '/genome.fasta')
    if len(fasta_file) != 1:
        fasta_file_gz = glob.glob(REFGENOMEFOLDER + '/genome.fasta.gz')
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOMEFOLDER)
        else:
            refgenome=fasta_file_gz[0]
    else:
        refgenome=fasta_file[0]
    refnt = np.zeros(position.shape[0],dtype=object)
    pos_counter=0
    for index,record in enumerate(SeqIO.parse(refgenome, format='fasta')):
        poschr = position[ position[:,0]==index+1 , 1]
        for sglpos in poschr:
            refnt[pos_counter]=record.seq[int(sglpos)]
            pos_counter += 1
    return refnt
    
def ana_mutation_quality(Calls,Quals):
    # This functions aims at providing the FQ value for every SNP position 
    # Across all pairwise different allele calls, it reports the best FQ value among the minimum FQ values per pair
    # NOTE: This function requires some more efficient coding!
    [Nmuts, NStrain] = Calls.shape ;
    MutQual = np.zeros((Nmuts,1)) ; 
    MutQualIsolates = np.zeros((Nmuts,2)); 
    
    # generate template index array to sort out strains gave rise to reported FQ values
    s_template=np.zeros( (len(Calls[0,:]),len(Calls[0,:])) ,dtype=object)
    for i in range(s_template.shape[0]):
        for j in range(s_template.shape[1]):
            s_template[i,j] = str(i)+"_"+str(j)

    for k in range(Nmuts):
        if len(np.unique(np.append(Calls[k,:], 4))) <= 2: # if there is only one type of non-N (4) call, skip this location
            MutQual[k] = np.nan ;
            MutQualIsolates[k,:] = 0; 
        else:
            c = Calls[k,:] ; c1 = np.tile(c,(c.shape[0],1)); c2 = c1.transpose() # extract all alleles for pos k and build 2d matrix and a transposed version to make pairwise comparison
            q = Quals[k,:] ; q1 = np.tile(q,(q.shape[0],1)); q2 = q1.transpose() # -"-
            g = np.all((c1 != c2 , c1 != 4 , c2 != 4) ,axis=0 )  # no data ==4; boolean matrix identifying find pairs of samples where calls disagree (and are not N) at this position
            #positive_pos = find(g); # numpy has no find; only numpy where, which does not flatten 2d array that way
            # get MutQual + logical index for where this occurred
            MutQual[k] = np.max(np.minimum(q1[g],q2[g])) # np.max(np.minimum(q1[g],q2[g])) gives lower qual for each disagreeing pair of calls, we then find the best of these; NOTE: np.max > max value in array; np.maximum max element when comparing two arryas
            MutQualIndex = np.argmax(np.minimum(q1[g],q2[g])) # return index of first encountered maximum!
            # get strain ID of reorted pair (sample number)
            s = s_template
            strainPairIdx = s[g][MutQualIndex]
            MutQualIsolates[k,:] = [strainPairIdx.split("_")[0], strainPairIdx.split("_")[1]]
            
    return [MutQual,MutQualIsolates]

def major_allele(arr):
    ''' returns 1-dimensional array of size arr.shape[0] with the major allele index [0:3] or NA [4] (if ambigous)'''
    # NOTE: if major NT ambigous (multiple alleles with same number of occurence I report allele with lower numeric value. could be improved as it might lead to negligible (imo) bias). Also loop could prob be removed.
    # input 2d arr (2nd dim can be 1!) with nucleotide indices (0:4)
    nonNA_out = (arr != 4)
    out_idx = []
    for i,row in enumerate(arr):
        if np.any(nonNA_out[i,:]): # any majorNT found
            row = row[nonNA_out[i,:]]
            row_ct = np.unique(row,return_counts=True)
            idx = np.where(row_ct[1] == np.max(row_ct[1]) )[0]
            out_idx.append(row_ct[0][idx][0]) # if multiple outgroup alleles same count I take the first. Could be changed to NA or add more outgroup samples for refined inference.
        else: # no majorNT
            out_idx.append(4)
    out_idx = np.array(out_idx)
    print(np.sum(out_idx == 4),'/', out_idx.size ,'elements of p have no major allele (ie. 4)!'  )
    return out_idx


def div_major_allele_freq(cnts,return_all_minorAF=False):
    # define matrices with major allele freq; major Nucl. (maNT); minor Nucl.; minor AF
    
    c=cnts[:,0:4,:]+cnts[:,4:8,:]; # flatten frw and rev ATCG counts    

    sorted_arr = np.sort(c,axis=1) # return sorted_arr matrix; axis 0 in 3d array == sort by col
    sortedpositions = np.argsort(c,axis=1) # return matrix indices of sort;axis 0 in 3d array == sort by col
    
    maxcount = sorted_arr[:,3:4:,:] # get allele counts for major allele (4th row); weird "3:4:" indexing required to maintain 3d structure
    minorcount = sorted_arr[:,2:3:,:] # get allele counts for first minor allele (3rd row); tri/quadro-allelic ignored; weird "2:3:" indexing required to maintain 3d structure and get 3rd row
    all_minorcount = sorted_arr[:,0:3:,:]

    with np.errstate(divide='ignore', invalid='ignore'): # suppress warning for division by 0
        maf = maxcount / sorted_arr.sum(axis=1,keepdims=True)
        minorAF = minorcount / sorted_arr.sum(axis=1,keepdims=True)
        all_minorAF = all_minorcount / sorted_arr.sum(axis=1,keepdims=True)

    maf = np.squeeze(maf,axis=1) # turn 2d; axis=1 to keep 2d structure when only one position!
    maf[ np.isnan(maf) ] = 0 # set to 0 to indicate no data
    minorAF = np.squeeze(minorAF,axis=1) 
    minorAF[ np.isnan(minorAF) ] = 0 # set to 0 to indicate no data/no minor AF
    all_minorAF = np.squeeze(all_minorAF)
    all_minorAF[ np.isnan(all_minorAF) ] = 0
    
    majorNT = np.squeeze(sortedpositions[:,3:4:,:],axis=1) # index position in sortedpositions represents allele position ATCG; axis=1 to keep 2d structure when only one position!
    minorNT = np.squeeze(sortedpositions[:,2:3:,:],axis=1) # -"-

    # Note: If counts for all bases are zero, then sort won't change the order
    # (since there is nothing to sort), thus majorNT.minorNT put to 4 (ie NA) using maf (REMEMBER: minorAF==0 is a value!)

    majorNT[maf==0] = 4
    minorNT[maf==0] = 4    
    
    if return_all_minorAF:
        return [maf.transpose(), majorNT.transpose(), minorNT.transpose(), minorAF.transpose(), all_minorAF.transpose()]

    return [maf.transpose(), majorNT.transpose(), minorNT.transpose(), minorAF.transpose()]



###################################################################################################
## # # # # # # # # # # # # # # # # All Filtering Functions # # # # # # # # # # # # # # # # # # # ##
###################################################################################################
# IO
def parse_bed_zero_covg_regions(path_to_bed_zero_covg_covg, p, scafNames, chrStarts, cutoff):
    """
    Generates an array of positions that should be masked due to falling within coverage islands.
    Reach to get to the next coverage island defined elsewhere (snakemake, other script) but cutoff
    defines % of uncovered region which must be uncovered to mask any covered positions within.

    e.g., covg = =======
         noncovered = _
         genome: ====___________________=_____________

         remove the isolated = covered position as likely contamination
    """
    # Use gzip.open for gzipped files
    print(f'Processing bed zero coverage file: {path_to_bed_zero_covg_covg}')
    with gzip.open(path_to_bed_zero_covg_covg, mode='rt') as file:  # Note the 'rt' mode for reading text
        lines = csv.reader(file, delimiter="\t")
        output_set = set()
        for line in lines:
            if float(line[6]) > cutoff:
                to_add_to_range = chrStarts[np.where(scafNames == line[0])[0]][0]
                output_set.update([x for x in range(int(line[1]) + to_add_to_range, int(line[2]) + 1 + to_add_to_range)])
    should_be_masked = []
    for index, pos in enumerate(p):
        if pos in output_set:
            should_be_masked.append(index)
    should_be_masked = np.array(should_be_masked)
    return should_be_masked

# filtering
def create_ranges(test_p,dist_to_run):
    # much more efficient version of old check for recombination range creation
    # including automatic ignoring of redundant ranges
    # runs in ~O(N) time
    current_index_limit_leading=0
    current_index_limit_trailing=0
    range_to_check_set=set()
    max_p_index=len(test_p)
    for position in test_p:
        # iterate up trailing limit
        lower_lim=position-dist_to_run
        while test_p[current_index_limit_trailing] < lower_lim:
            current_index_limit_trailing +=1
        # iterate up leading limit
        upper_lim=position+dist_to_run
        while test_p[current_index_limit_leading] < upper_lim:
            if current_index_limit_leading == max_p_index-1:
                break
            else:
                current_index_limit_leading += 1
        range_to_check_set.add((current_index_limit_trailing,current_index_limit_leading))
    return range_to_check_set


def findrecombinantSNPs(p, mutantAF, distance_for_nonsnp, corr_threshold_recombination,call_matrix_failing_qc=[]):
    # find ranges:
    ranges_to_search=create_ranges(p,distance_for_nonsnp)
    # prepare output and check that ranges exist
    nonsnp_idx = np.zeros(0,dtype='int')
    if len(ranges_to_search) == 0:
        return nonsnp_idx,np.zeros(p.shape)
    # iterate over range to ID highly correlated SNPs
    for range_start,range_end in ranges_to_search:
        r = mutantAF[[range_start,range_end],:]
        if len(call_matrix_failing_qc) > 0:
            calls_masked_any_sample=(np.sum(call_matrix_failing_qc[[range_start,range_end],:],axis=0)>0)
            r_only_unmasked_samples=r[:,~calls_masked_any_sample]
            corrmatrix = np.corrcoef(r_only_unmasked_samples)
        else:
            corrmatrix = np.corrcoef(r)
        [a,b]=np.where(corrmatrix > corr_threshold_recombination)
        nonsnp_idx=np.concatenate((nonsnp_idx,np.arange(range_start,range_end)[a[np.where(a!=b)]]))
    nonsnp_idx=np.unique(nonsnp_idx)
    ## print number of sites which have indication for recombination
    print('\n' + str(nonsnp_idx.shape[0]) + ' of a total ' + str(p.shape[0]) + ' ('  + str((len(nonsnp_idx)/len(p))*100) + '%) positions in goodpos were found to be recombinant.')
    nonsnp_bool = np.zeros(p.shape) # build bool array of length p with nonsnps == True
    if nonsnp_idx.size > 0:
        nonsnp_bool[nonsnp_idx] = 1
    nonsnp_bool = nonsnp_bool.astype(bool)
    return [nonsnp_idx, nonsnp_bool]

def filter_bed_cov_hist(bed_path,p,scafNames,chrStarts,sampleNames,coverage,cutoff,two_tailed=False,upper=True):
    """
    NOTE: For ancient DNA quality control!!!
    outputs a boolean matrix of index of [p]x[samplenames] to mask from basecall, due to being in the top coverage percentile when using bedtools coverage histogram (or bottom if upper=False, or if two_tailed=True)"""
    # get paths
    bed_histogram_files = glob.glob(bed_path)
    # index accounting
    p_chr_indices=[0]+[np.max(np.where(p < x + 1))+1 for x in chrStarts[1:]]
    # generate output matrix
    to_be_masked_array_covg_percentile=np.full(coverage.shape,False)
    for bed_index,bed_hist in enumerate(bed_histogram_files):
        this_sample_index=-1
        print(f'Processing bed covg hist file: {bed_hist}')
        for index,samplename in enumerate(sampleNames):
            if samplename in bed_hist:
                this_sample_index=index
        if this_sample_index==-1:
            print(f'Warning: Bedfile {bed_hist} does not have a corresponding samplename in sampleNames input array, skipping.')
        else:
            this_sample_name_cutoffs=cutoff_bed_covg_histograms(bed_hist,cutoff,two_tailed,upper)
            print(this_sample_name_cutoffs)
            ## now, mask positions above the percentile, by individual chrom cutoffs
            for index,scaf in enumerate(scafNames):
                this_sample_name_cutoffs_this_scaf=this_sample_name_cutoffs[scaf][1]
                if index<len(scafNames)-1:
                    start=p_chr_indices[index]
                    end=p_chr_indices[index+1]
                    p_to_include_this_chrom=np.array(range(start,end)) ## true/false 
                else:
                    p_to_include_this_chrom=np.array(range(start,len(p)))
                to_mask=p_to_include_this_chrom[np.in1d(p_to_include_this_chrom,np.where(coverage[:,this_sample_index] > this_sample_name_cutoffs_this_scaf)[0])]
                to_be_masked_array_covg_percentile[to_mask,this_sample_index]=True
    return to_be_masked_array_covg_percentile

def filter_bed_0_cov_regions(bed_path,p,scafNames,chrStarts,sampleNames,cutoff):
    """outputs a boolean matrix of index of [p]x[samplenames] to mask from basecall, due to being in regions of 0 coverage (eg are coverage islands)"""
    # get paths
    bed_zero_covg_files=glob.glob(bed_path)
    # generate output matrix
    to_be_masked_array_0_covg_regions=np.full((len(p),len(sampleNames)),False)
    for bed_index,bed_zero in enumerate(bed_zero_covg_files):
        this_sample_index=-1
        for index,samplename in enumerate(sampleNames):
            if samplename in bed_zero:
                this_sample_index=index
        if this_sample_index==-1:
            print(f'Warning: Bedfile {bed_zero} does not have a corresponding samplename in sampleNames input array, skipping.')
        else:
            to_mask=parse_bed_zero_covg_regions(bed_zero,p,scafNames,chrStarts,cutoff)
            ## now, mask positions above the percentile 
            to_be_masked_array_0_covg_regions[to_mask,this_sample_index]=True
    return to_be_masked_array_0_covg_regions

def cutoff_bed_covg_histograms(path_to_bed_covg_hist, cutoff, two_tailed=True, upper=True):
    """generates dictionary of scafnames --> cutoff values for coverage, based on bedtools covg histogram"""
    # get cutoffs for coverage
    if cutoff >= 0.5: cutoff = 1 - cutoff
    
    # Check if the file is gzipped based on its extension
    if path_to_bed_covg_hist.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    
    with opener(path_to_bed_covg_hist, 'rt') as file:  # 'rt' mode for reading text from a gzipped file
        lines = csv.reader(file, delimiter="\t")
        output_hist = {}
        for line in lines:
            if line[0] not in output_hist:
                output_hist[line[0]] = {line[1]: line[4]}
            else:
                output_hist[line[0]][line[1]] = line[4]
    
    if two_tailed:
        cutoffs = (cutoff, 1 - cutoff)
    else:
        if upper:
            cutoffs = (0, 1 - cutoff)
        else:
            cutoffs = (cutoff, 1)
    
    cutoffs_scafs = {scaf: [0, 0] for scaf in output_hist.keys()}
    for scaf in output_hist:
        cumulative_perctile = 0
        for cov, value in output_hist[scaf].items():
            cov, value = float(cov), float(value)
            cumulative_perctile += value
            if cumulative_perctile < cutoffs[0]:
                cutoffs_scafs[scaf][0] = cov
            if cumulative_perctile > cutoffs[1]:
                cutoffs_scafs[scaf][1] = cov
                break
    return cutoffs_scafs

def find_calls_near_heterozygous_sites(p, minorAF, distance_to_check, heterozygosity_threshold, indices_of_ancient_samples=np.array([])):
    """
    Retruns highly correlated positoins, with optional masking of positions that fail some QC check (eg not enough covg --> call of N)
    NOTE: By default, SNPs called N will be factored into the correlation matrix, not recommended with ancient data (!!)
    hackathon add default behavior (either by using calls matrix (but does covar matrix work with that? or just has mut and calls matrix))
    """
    #look for recombination regions
    if len(indices_of_ancient_samples) == 0:
        print('Conducting nearby heterozygosity check on all samples, note: this filter should only be applied to ancient samples!')
        indices_of_ancient_samples=np.array([x for x in range(0,minorAF.shape[1])])
    site_samples_to_mask = np.full(minorAF.shape,False)
    ranges_to_search=create_ranges(p,distance_to_check)
    for range_start,range_end in ranges_to_search:
        r = minorAF[[range_start,range_end],:]
        for sample_index in indices_of_ancient_samples:
            if np.any(r[:,sample_index]>heterozygosity_threshold):
                site_samples_to_mask[[range_start,range_end],sample_index]=True
    return site_samples_to_mask

