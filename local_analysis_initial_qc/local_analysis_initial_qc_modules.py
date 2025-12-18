
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

# Processing of gff
"""def parse_gff(REFGENOMEFOLDER,ScafNames,ortholog_info_series=pd.Series(dtype = 'int64'),forceReDo=False):
    # parse gff file (tested with version3) with genome annotation (gff or gff.gz)
    # requires one gff file in REFGENOMEFOLDER, which is detected automatically
    # provide ScafNames to maintain same order of contigs/chr/scaffolds in output list as in previously generated ref-genome-based variables, eg. contig_positions,ChrStarts, GenomeLength
    # some columns potentially contain multiple entries. Those are separated by ";"
    # no data is always reported as "."
    # https://biopython.org/wiki/GFF_Parsing 
    # NOTE: annotation is read from gff file!
    # only execute if dataframes not yet made
    if os.path.isfile(REFGENOMEFOLDER+"/annotation_genes.pandas.py.pk1") and (forceReDo == False):
        afile = open(REFGENOMEFOLDER+"/annotation_genes.pandas.py.pk1", 'rb')
        list_of_dataframes = pickle.load(afile)
        afile.close()
        return list_of_dataframes
    else: # search for gff or gff.gz
        examiner = GFF.GFFExaminer()        
        # remake examiner which just checks what is inside the gff
        gff_file = glob.glob(REFGENOMEFOLDER + '/*.gff*') # search for gff or gff.gz or gff3
        if len(gff_file) != 1:
            raise ValueError('Either no gff(.gz) file or more than 1 *gff(.gz) file found in ' + REFGENOMEFOLDER)
        if gff_file[0][-2:] == 'gz':
            encoding = 'gzip'
        else:
            encoding = 'unzip'
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open # define opening procedure gzip.open or open
        with _open(gff_file[0]) as gff_handle:
            possible_limits = examiner.available_limits(gff_handle)
        # start processing
        chromosomes = possible_limits["gff_id"].keys() # list of tuples containing contig ID; do not loop over those, as I only loop over contigs observed (ScafNames)...might need a fix if !=
        tagnumber_counter = 0 # unique, continous numerical id for all features across all chromosomes
        list_of_dataframes = [] # output: each element contains a dtaframe w/ all annotation info for chr. ordered as in ScafNames.
        # define gff_type to extracted. Use all but gene and region. gene annotated in NCBI but lacks most info that extra CDS has. Region just describes chromosome extensions
        limits = dict(gff_type = [i[0] for i in possible_limits['gff_type'].keys() if i[0] != 'gene' and i[0] != 'region'] )
        # loop over chr with established order
        for chrom in ScafNames:
            limits["gff_id"] = [chrom]
            # limits = dict(gff_id=[chrom])
            with _open(gff_file[0]) as gff_handle:
                for rec in GFF.parse(gff_handle, limit_info=limits):
                    # for loop over every chr but only defined [limits] has data. Memory efficient!
                    if rec.id == chrom:
                        # if chr has any feature build list of dicts and append to list_of_dataframes, else append empty dataframe
                        if len(rec.features) > 0:
                            # test if seq object part of gff (prokka-based yes, but NCBI-based no >> then load ref genome.fasta)
                            if len(rec.seq) == rec.seq.count('?'):
                                for seq_record in SeqIO.parse(REFGENOMEFOLDER+"/genome.fasta", "fasta"):
                                    if seq_record.id == rec.id:
                                        rec.seq = seq_record.seq
                                if len(rec.seq) == rec.seq.count('?'): # test if succesful
                                    print('Warning: No reference genome found that matches chromosome:' + rec.id)
                            print(rec.id)
                            lod_genes = [] # list-of-dictionary. Easy to turn to pd.dataframe
                            for gene_feature in rec.features:
                                gene_dict = {}
                                tagnumber_counter += 1
                                gene_dict['type'] = gene_feature.type
                                gene_dict['locustag'] = gene_feature.id
                                # add ortholog info if locustag (eg. repeat region has none)
                                #print("gene_feature.id   "+gene_feature.id)
                                if gene_feature.id != "" and gene_feature.type == 'CDS' and not ortholog_info_series.empty:
                                    gene_dict['orthologtag'] = ortholog_info_series[ortholog_info_series.str.findall(gene_feature.id).str.len() == 1].index[0]
                                #print(rec.id+"   "+gene_dict['locustag'])
                                if 'gene' in gene_feature.qualifiers.keys():
                                    gene_dict['gene'] = ";".join(gene_feature.qualifiers['gene'])
                                else:
                                    gene_dict['gene'] = "." # add "." instead of []
                                gene_dict['type'] = gene_feature.type
                                gene_dict['locustag'] = gene_feature.id
                                if gene_dict['type'] == "CDS" or gene_dict['type'] == "gene":
                                    gene_dict['tagnumber'] = tagnumber_counter
                                else:
                                    gene_dict['tagnumber'] = 0
                                if 'product' in gene_feature.qualifiers.keys():
                                    gene_dict['product'] = ";".join(gene_feature.qualifiers['product']) 
                                else:
                                    gene_dict['product'] = "."
                                if 'protein_id' in gene_feature.qualifiers.keys():
                                    gene_dict['protein_id'] = gene_feature.qualifiers['protein_id']
                                else:
                                    gene_dict['protein_id'] = "."
                                if "Dbxref" in gene_feature.qualifiers.keys():
                                    gene_dict['db_xref'] = ";".join(gene_feature.qualifiers['Dbxref'])
                                else:
                                    gene_dict['db_xref'] = "."
                                if "note" in gene_feature.qualifiers.keys():
                                    gene_dict['note'] = ";".join(gene_feature.qualifiers['note'])
                                elif "Note" in gene_feature.qualifiers.keys():
                                    gene_dict['note'] = ";".join(gene_feature.qualifiers['Note'])
                                else:
                                    gene_dict['note'] = "."
                                if 'phase' in gene_feature.qualifiers.keys():
                                    gene_dict['codon_start'] = int(gene_feature.qualifiers['phase'][0])
                                else:
                                    gene_dict['codon_start'] = "."
                                gene_dict['indices'] = [gene_feature.location.start.position,gene_feature.location.end.position]
                                gene_dict['loc1'] = gene_feature.location.start.position # automatically 0-based
                                gene_dict['loc2'] = gene_feature.location.end.position
                                gene_dict['location'] = gene_feature.location
                                gene_dict['strand'] = gene_feature.location.strand
                                dna_seq = rec.seq[gene_feature.location.start:gene_feature.location.end]
                                if gene_dict['strand'] == 1:
                                    gene_dict['sequence'] = dna_seq
                                elif gene_dict['strand'] == -1:
                                    gene_dict['sequence'] = dna_seq.reverse_complement()
                                else:
                                    gene_dict['sequence'] = dna_seq # eg. repeat region
                                # translation, add info where codon starts if info was available. Usually it starts at 0
                                if isinstance( gene_dict['codon_start'] , int):
                                    sequence2translate = gene_dict['sequence'][gene_dict['codon_start']:]
                                    gene_dict['translation'] = sequence2translate.translate(table="Bacterial") # bacterial genetic code GTG is a valid start codon, and while it does normally encode Valine, if used as a start codon it should be translated as methionine. http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:translation
                                elif gene_dict['type'] == "CDS":
                                    sequence2translate = gene_dict['sequence']
                                    gene_dict['translation'] = sequence2translate.translate(table="Bacterial")
                                else:
                                    gene_dict['translation'] = "." # all non-CDS (RNA's or repeat regions) not translated (as those are sometimes also off-frame)
                                lod_genes.append(gene_dict)
                            # after all features parsed > add as dataframe to output list
                            # ATTENTION: SORT pandas DF (annotation not necessarily sorted)
                            df_sort = pd.DataFrame(lod_genes)
                            df_sort = df_sort.sort_values(by=['loc1'])
                            list_of_dataframes.append(df_sort)
                            
                        else:
                            list_of_dataframes.append( pd.DataFrame([{}]) )
        afile = open(REFGENOMEFOLDER+"/annotation_genes.pandas.py.pk1", 'wb')
        pickle.dump(list_of_dataframes, afile)
        afile.close()
        return list_of_dataframes
"""

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

###################################################################################################
## # # # # # # # # # # # # # All Annotation Mutation Functions # # # # # # # # # # # # # # # # # ##
###################################################################################################

def annotate_mutations(annotation_genes , p_gp , refnti_gp , ancnti_gp , calls_gp , counts_gp , hasmutation_gp , mutQual_gp, promotersize, ref_genome_folder, goodpos2use=[]):
    ''' produce dataframe with annotation info for each goodposition used for tree
    all positions are 1-based! (pos, nt_pos, loc1/2, distance1/2)
    # function combines TDL: annotate_mutations_gb.m and append_annotations.m 
    CHANGELOG:
    17.10.2022 IL: Added parsing of ancestral AA state by getting annotation info from ancnti_gp'''
    # p_gp: genomic position of goodpos
    # p_gp=p[goodpos2use]
    # refnti_gp=refnti_m[goodpos2use,:]
    # calls_gp=calls[goodpos2use,:]
    # counts_gp=counts[:,:,goodpos2use]
    # hasmutation_gp=hasmutation[goodpos2use,:]
    # mutQual_gp = mutQual[goodpos2use,].flatten() 
    # ref_genome_folder=REFGENOMEFOLDER
    if len(goodpos2use) > 0:
        p_gp = p_gp[goodpos2use]
        refnti_gp = refnti_gp[goodpos2use,:]
        ancnti_gp = ancnti_gp[goodpos2use,:]
        calls_gp = calls_gp[goodpos2use,:]
        counts_gp = counts_gp[:,:,goodpos2use]
        hasmutation_gp = hasmutation_gp[goodpos2use,:]
        mutQual_gp = mutQual_gp[goodpos2use,].flatten()

    [maf_gp, maNT_gp, minorNT_gp, minorAF_gp] = div_major_allele_freq(counts_gp)
    nts = ['A','T','C','G'] #'atcg';
    rc = ['T','A','G','C']  # 'tagc';
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    [locus_tag,chr_tag] = genomic_position_all(annotation_genes, genomelength, chrstarts); #loc_tag numbers all genes across genome (intergenic:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    mut_genenum = chr_tag[p_gp]
    loc_tag_p_gp = locus_tag[p_gp]
    pos_chr = p2chrpos(p_gp,chrstarts) # get two-col table chr,pos   
    lod_mutAnno = []
    for i,pos in enumerate(p_gp):
        #print(i,pos)
        mut_annotations = {}
        mut_annotations['gene_num'] = mut_genenum[i]
        mut_annotations['gene_num_global'] = loc_tag_p_gp[i]
        mut_annotations['chr'] = pos_chr[i][0]
        mut_annotations['pos'] = pos_chr[i][1] + 1 # turn 0-based pos into 1-based pos.
        mut_annotations['p'] = p_gp[i]
        mut_annotations['quals'] = mutQual_gp[i]
        mut_annotations['nts'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_ref'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_anc'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_alt'] = "." # default. overwritten below if NT defined.
        p_chr = pos_chr[i][0]
        if mut_genenum[i] == int(mut_genenum[i]): # intragenic
            p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])-1] # both -1 necessary bcs list of df and df 0-based
            mut_annotations['product'] = p_anno.loc['product']
            mut_annotations['gene'] = p_anno.loc['gene']
            mut_annotations['protein_id'] = p_anno.loc['protein_id']
            mut_annotations['strand'] = p_anno.loc['strand']
            mut_annotations['loc1'] = p_anno.loc['loc1'] + 1 # 1-based first position
            mut_annotations['loc2'] = p_anno.loc['loc2'] # last position of gene (inclusive)
            mut_annotations['sequence'] = p_anno.loc['sequence']
            mut_annotations['protein_id'] = p_anno.loc['protein_id']
            mut_annotations['note'] = p_anno.loc['note']
            mut_annotations['locustag'] = p_anno.loc['locustag']
            if 'orthologtag' in p_anno:
                mut_annotations['orthologtag'] = p_anno.loc['orthologtag']
            mut_annotations['translation'] = p_anno.loc['translation']
            if p_anno.loc['strand'] == 1: # get position within gene, consider strandedness
                mut_annotations['nt_pos'] = mut_annotations['pos'] - mut_annotations['loc1']+1 # pos/loc1 1-based. nt_pos 1-based. +1 to get 1-based nt_pos in gene (checked)!
            elif p_anno.loc['strand'] == -1:
                mut_annotations['nt_pos'] = p_anno.loc['loc2'] - mut_annotations['pos'] +1 # pos/loc2 1-based. +1 to get 1-based nt_pos in gene (checked)!
            else:
                mut_annotations['nt_pos'] = "." # can happen. eg. Crispr
            if mut_annotations['nt_pos'] == ".": # Observed with special 'type's like crispr annotated. rare! leads to no AA inference.
                aan = 9999999
            else:
                aan = int( (mut_annotations['nt_pos']-1 )/3 ) + 1 # get #codon that harbours mutation. 1-based.
                ncn = mut_annotations['nt_pos'] - ((aan-1)*3) # get #nucl within triplett. 1-based
                mut_annotations['aa_pos'] = aan
            codons_ls = []
            aa_ls = [] #potential AA changes at given codon
            if len(mut_annotations['sequence']) >= (aan*3) and mut_annotations['translation'] != ".":
                codon = mut_annotations['sequence'][aan*3-2-1:aan*3] # -1 bcs seq-object 0-based; but positional argument aan not 
                codon = [n for n in codon  ] # turn seq.object to list, seq object does not allow reassignment
                for n in range(4): # test all four nucleotide options for resulting AA change
                    if p_anno.loc['strand'] == 1:
                        codon[ncn-1] = nts[n]
                    else:
                        codon[ncn-1] = rc[n]
                    codonSeqO = Seq( "".join(codon)) 
                    codons_ls.append(codonSeqO)
                    aa_ls.append(codonSeqO.translate())
            mut_annotations['codons'] = codons_ls
            mut_annotations['AA'] = [aa[0] for aa in aa_ls]
            # append_annotations intragenic
            if len(mut_annotations['AA']) < 4:
                mut_annotations['type'] = 'U'
            mut_annotations['AA_gt'] = '' # Fill in annotations with whether NS or Syn mutation
            mut_annotations['AA_gt_anc'] = '' 
            if np.unique(refnti_gp[i,:])[0] != 4: # only if ref defined; [0] ok...see next line comment
                mut_annotations['nt_ref'] = nts[np.unique(refnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                mut_annotations['nts'] = mut_annotations['nt_ref']
                if np.unique(ancnti_gp[i,:])[0] != 4:
                    mut_annotations['nt_anc'] = nts[np.unique(ancnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                if len(mut_annotations['AA']) == 4:
                    mut_annotations['AA_gt'] = mut_annotations['AA_gt'] + mut_annotations['AA'][ np.unique(refnti_gp[i,:])[0] ] # the AA list order corresponds to NTS list!
            # extract derived genotype(s) and according AA across all samples
            # call Nonsynonymous if >1 amino acid type exists across all called samples (eg ref AA and alt )
            for j,callidx in enumerate(calls_gp[i,:]):
                #print(str(j)+" "+str(callidx))
                #fixed mutation
                # i is index of current pos, j is index of samples 
                if hasmutation_gp[i,j] == True :
                    if calls_gp[i,j] < 4:
                        # add each NT called, later only keep uniq, but kept order
                        mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_gp[i,j]]
                        mut_annotations['nt_alt'] = nts[calls_gp[i,j]]
                        if len(mut_annotations['AA']) == 4:
                            mut_annotations['AA_gt'] = mut_annotations['AA_gt'] + mut_annotations['AA'][ maNT_gp[i,j] ]
                            if np.unique(ancnti_gp[i,:])[0] != 4:
                                mut_annotations['AA_gt_anc'] = mut_annotations['AA_gt_anc'] + mut_annotations['AA'][ np.unique(ancnti_gp[i,:])[0] ]
            if len(mut_annotations['AA']) == 4:
                mut_annotations['type'] = 'S' # eventually overwritten below if N
            # remove duplicates
            mut_annotations['AA_gt_anc'] = ''.join(OrderedDict.fromkeys( mut_annotations['AA_gt_anc'] ).keys()) # get only unique AA and keep order
            mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
            mut_annotations['AA_gt'] = ''.join(OrderedDict.fromkeys( mut_annotations['AA_gt'] ).keys()) # get only unique AA and keep order
            # record if nonsynonymous mutation
            if len(mut_annotations['AA_gt'])>1:
                mut_annotations['type'] = 'N'
            # Record all observed mutations across all isolates; E.g. K134Y, W47A, etc.
            if len(mut_annotations['AA_gt'])>1:
                mut_annotations['muts'] = []
                ancAA = mut_annotations['AA_gt'][0]
                derAAs = mut_annotations['AA_gt'][1:]
                for j,derAA in enumerate(derAAs):
                    mut_annotations['muts'].append( ancAA+str(mut_annotations['aa_pos'])+derAA )
            else:
                mut_annotations['muts'] = "."
            mut_annotations['NonSyn_rel_ref'] = len(mut_annotations['AA_gt'])>1
            if len(mut_annotations['AA_gt_anc'])>1:
                mut_annotations['muts_anc'] = []
                ancestralAA = mut_annotations['AA_gt_anc'][0]
                derivedAAs = mut_annotations['AA_gt_anc'][1:]
                for j,derAA in enumerate(derivedAAs):
                    mut_annotations['muts_anc'].append( ancestralAA+str(mut_annotations['aa_pos'])+derAA )
            else:
                mut_annotations['muts_anc'] = "."
            if mut_annotations['AA_gt'] != mut_annotations['AA_gt_anc']:
                mut_annotations['NonSyn_rel_anc'] = True
            else: 
                mut_annotations['NonSyn_rel_anc'] = False
            if mut_annotations['NonSyn_rel_anc'] or mut_annotations['NonSyn_rel_ref']:
                mut_annotations['NonSyn'] = True
            else:
                mut_annotations['NonSyn'] = False
        else: #intergenic
            if int(mut_genenum[i])>0: # get info for gene prior SNP (if any)
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])-1] # both -1 necessary bcs list of df and df 0-based
                mut_annotations['gene1'] = p_anno.loc['gene']
                mut_annotations['locustag1'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag1'] = p_anno.loc['orthologtag']
                mut_annotations['product1'] = p_anno.loc['product']
                mut_annotations['distance1'] = mut_annotations['pos'] - p_anno.loc['loc2']
                if p_anno.loc['strand'] == -1:
                    mut_annotations['distance1'] = mut_annotations['distance1'] * -1

            if int(mut_genenum[i]+0.5) <= annotation_genes[p_chr-1].shape[0] and annotation_genes[p_chr-1].shape[1] !=0: # get info gene after SNP (if any); second conditional to evade empty chr
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])] # -1 necessary bcs list of df 0-based; gene_id 0-based by we want following
                mut_annotations['gene2'] = p_anno.loc['gene']
                mut_annotations['locustag2'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag2'] = p_anno.loc['orthologtag']
                mut_annotations['product2'] = p_anno.loc['product']
                mut_annotations['distance2'] = p_anno.loc['loc1'] - mut_annotations['pos'] +1 # +1 to get correct bp distance
                if p_anno.loc['strand'] == 1:
                    mut_annotations['distance2'] = mut_annotations['distance2'] * -1
            # append_annotations intragenic
            #print(mut_annotations['distance1'])
            if ( 'distance1' in mut_annotations and mut_annotations['distance1'] > (-1*promotersize) and mut_annotations['distance1'] < 0) or ( 'distance2' in mut_annotations and mut_annotations['distance2'] > (-1*promotersize) and mut_annotations['distance2'] < 0):
                mut_annotations['type'] = 'P'
            else:
                mut_annotations['type'] = 'I'
            # define nts (repeat of intragenic code)
            if np.unique(refnti_gp[i,:])[0] != 4: # only if ref defined; [0] ok...see next line comment
                mut_annotations['nt_ref'] = nts[np.unique(refnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                mut_annotations['nts'] = mut_annotations['nt_ref']
                if np.unique(ancnti_gp[i,:])[0] != 4:
                    mut_annotations['nt_anc'] = nts[np.unique(ancnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
            # extract derived genotype(s) across all samples
            for j,callidx in enumerate(calls_gp[i,:]):
                #print(str(j)+" "+str(callidx))
                #fixed mutation
                if hasmutation_gp[i,j] == True:
                    if calls_gp[i,j] < 4:
                        mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_gp[i,j]]
                        mut_annotations['nt_alt'] = nts[calls_gp[i,j]]
            # remove duplicates
            mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
        lod_mutAnno.append(mut_annotations)
    dataframe_mut = pd.DataFrame(lod_mutAnno)
    return dataframe_mut

def annotate_mutations_single_sample_vs_custom_nt_call(annotation_genes , p_goodpos , custom_nt_calls_goodpos , calls_single_sample_goodpos , hasmutation_goodpos , promotersize, ref_genome_folder, goodpos=[]):
    ''' 
    produces annotation_mutations-like dataframe with annotation info for each goodposition used for tree for a given sample
    Comparison and type called based on "calls" from user provided goodpos
    all positions are 1-based! (pos, nt_pos, loc1/2, distance1/2)
    CHANGELOG:
    17.10.2022 IL: Added parsing of ancestral AA state by getting annotation info from ancnti_gp
    23.03.2023 IL: Removed reference nt submission for clarity, all reliance on provided custom_nt_call call'''
    # p_gp: genomic position of goodpos
    # p_gp=p[goodpos2use]
    # refnti_gp=refnti_m[goodpos2use,:]
    # calls_gp=calls[goodpos2use,:]
    # counts_gp=counts[:,:,goodpos2use]
    # hasmutation_gp=hasmutation[goodpos2use,:]
    # mutQual_gp = mutQual[goodpos2use,].flatten() 
    # ref_genome_folder=REFGENOMEFOLDER
    if len(goodpos) > 0:
        p_goodpos = p_goodpos[goodpos]
        custom_nt_calls_goodpos = custom_nt_calls_goodpos[goodpos,:]
        calls_single_sample_goodpos = calls_single_sample_goodpos[goodpos,:]
        hasmutation_goodpos = hasmutation_goodpos[goodpos,:]
    nts = ['A','T','C','G'] #'atcg';
    rc = ['T','A','G','C']  # 'tagc';
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    [locus_tag,chr_tag] = genomic_position_all(annotation_genes, genomelength, chrstarts); #loc_tag numbers all genes across genome (intergenic:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    mut_genenum = chr_tag[p_goodpos]
    loc_tag_p_gp = locus_tag[p_goodpos]
    pos_chr = p2chrpos(p_goodpos,chrstarts) # get two-col table chr,pos   
    lod_mutAnno = []
    for i,pos in enumerate(p_goodpos):
        #print(i,pos)
        mut_annotations = {}
        mut_annotations['gene_num'] = mut_genenum[i]
        mut_annotations['gene_num_global'] = loc_tag_p_gp[i]
        mut_annotations['chr'] = pos_chr[i][0]
        mut_annotations['pos'] = pos_chr[i][1] + 1 # turn 0-based pos into 1-based pos.
        mut_annotations['p'] = p_goodpos[i]
        mut_annotations['nts'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_ref'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_anc'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_alt'] = "." # default. overwritten below if NT defined.
        p_chr = pos_chr[i][0]
        if mut_genenum[i] == int(mut_genenum[i]): # intragenic
            p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])-1] # both -1 necessary bcs list of df and df 0-based
            mut_annotations['product'] = p_anno.loc['product']
            mut_annotations['gene'] = p_anno.loc['gene']
            mut_annotations['protein_id'] = p_anno.loc['protein_id']
            mut_annotations['strand'] = p_anno.loc['strand']
            mut_annotations['loc1'] = p_anno.loc['loc1'] + 1 # 1-based first position
            mut_annotations['loc2'] = p_anno.loc['loc2'] # last position of gene (inclusive)
            mut_annotations['sequence'] = p_anno.loc['sequence']
            mut_annotations['protein_id'] = p_anno.loc['protein_id']
            mut_annotations['note'] = p_anno.loc['note']
            mut_annotations['locustag'] = p_anno.loc['locustag']
            mut_annotations['chr_locustag'] = str(pos_chr[i][0]) + '_' + p_anno.loc['locustag']
            if 'orthologtag' in p_anno:
                mut_annotations['orthologtag'] = p_anno.loc['orthologtag']
            mut_annotations['translation'] = p_anno.loc['translation']
            if p_anno.loc['strand'] == 1: # get position within gene, consider strandedness
                mut_annotations['nt_pos'] = mut_annotations['pos'] - mut_annotations['loc1']+1 # pos/loc1 1-based. nt_pos 1-based. +1 to get 1-based nt_pos in gene (checked)!
            elif p_anno.loc['strand'] == -1:
                mut_annotations['nt_pos'] = p_anno.loc['loc2'] - mut_annotations['pos'] +1 # pos/loc2 1-based. +1 to get 1-based nt_pos in gene (checked)!
            else:
                mut_annotations['nt_pos'] = "." # can happen. eg. Crispr
            if mut_annotations['nt_pos'] == ".": # Observed with special 'type's like crispr annotated. rare! leads to no AA inference.
                aan = 9999999
            else:
                aan = int( (mut_annotations['nt_pos']-1 )/3 ) + 1 # get #codon that harbours mutation. 1-based.
                ncn = mut_annotations['nt_pos'] - ((aan-1)*3) # get #nucl within triplett. 1-based
                mut_annotations['aa_pos'] = aan
            codons_ls = []
            aa_ls = []
            if len(mut_annotations['sequence']) >= (aan*3) and mut_annotations['translation'] != ".":
                codon = mut_annotations['sequence'][aan*3-2-1:aan*3] # -1 bcs seq-object 0-based; but positional argument aan not 
                codon = [n for n in codon  ] # turn seq.object to list, seq object does not allow reassignment
                for n in range(4): # test all four nucleotide options for resulting AA change
                    if p_anno.loc['strand'] == 1:
                        codon[ncn-1] = nts[n]
                    else:
                        codon[ncn-1] = rc[n]
                    codonSeqO = Seq( "".join(codon)) 
                    codons_ls.append(codonSeqO)
                    aa_ls.append(codonSeqO.translate())
            mut_annotations['codons'] = codons_ls
            mut_annotations['AA'] = [aa[0] for aa in aa_ls]
            # append_annotations intragenic
            if len(mut_annotations['AA']) < 4:
                mut_annotations['type'] = 'U'
            mut_annotations['AA_gt_custom'] = '' 
            if custom_nt_calls_goodpos[i] != 4: # only if ref defined; [0] ok...see next line comment
                mut_annotations['nt_custom'] = nts[custom_nt_calls_goodpos[i]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                mut_annotations['nts'] = mut_annotations['nt_custom']
                if custom_nt_calls_goodpos[i] != 4:
                    mut_annotations['nt_custom'] = nts[custom_nt_calls_goodpos[i]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
            # extract derived genotype(s) and according AA across all samples
            # call Nonsynonymous if >1 amino acid type exists across all called samples (eg ref AA and alt )
            #print(str(j)+" "+str(callidx))
            #fixed mutation
            # i is index of current pos, j is index of samples 
            if hasmutation_goodpos[i] == True :
                if calls_single_sample_goodpos[i] < 4:
                    # add each NT called, later only keep uniq, but kept order
                    mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_single_sample_goodpos[i]]
                    mut_annotations['nt_alt'] = nts[calls_single_sample_goodpos[i]]
                    if len(mut_annotations['AA']) == 4:
                        mut_annotations['AA_gt_custom'] = mut_annotations['AA_gt_custom'] + mut_annotations['AA'][ calls_single_sample_goodpos[i] ]
                        if np.unique(custom_nt_calls_goodpos[i])[0] != 4:
                            mut_annotations['AA_gt_custom'] = mut_annotations['AA_gt_custom'] + mut_annotations['AA'][ custom_nt_calls_goodpos[i] ]
            if len(mut_annotations['AA']) == 4:
                mut_annotations['type'] = 'S' # eventually overwritten below if N
            # remove duplicates
            mut_annotations['AA_gt_custom'] = ''.join(OrderedDict.fromkeys( mut_annotations['AA_gt_custom'] ).keys()) # get only unique AA and keep order
            mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
            # record if nonsynonymous mutation
            if len(mut_annotations['AA_gt_custom'])>1:
                mut_annotations['type'] = 'N'
            # Record all observed mutations across all isolates; E.g. K134Y, W47A, etc.
            if len(mut_annotations['AA_gt_custom'])>1:
                mut_annotations['muts'] = []
                ancAA = mut_annotations['AA_gt_custom'][0]
                derAAs = mut_annotations['AA_gt_custom'][1:]
                for j,derAA in enumerate(derAAs):
                    mut_annotations['muts'].append( ancAA+str(mut_annotations['aa_pos'])+derAA )
            else:
                mut_annotations['muts'] = "."
            if mut_annotations['type'] == 'N':
                mut_annotations['NonSyn'] = True
            elif mut_annotations['type'] == 'S':
                mut_annotations['NonSyn'] = False
        else: #intergenic
            mut_annotations['NonSyn'] = np.NaN
            if int(mut_genenum[i])>0: # get info for gene prior SNP (if any)
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])-1] # both -1 necessary bcs list of df and df 0-based
                mut_annotations['gene1'] = p_anno.loc['gene']
                mut_annotations['locustag1'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag1'] = p_anno.loc['orthologtag']
                mut_annotations['product1'] = p_anno.loc['product']
                mut_annotations['distance1'] = mut_annotations['pos'] - p_anno.loc['loc2']
                if p_anno.loc['strand'] == -1:
                    mut_annotations['distance1'] = mut_annotations['distance1'] * -1

            if int(mut_genenum[i]+0.5) <= annotation_genes[p_chr-1].shape[0] and annotation_genes[p_chr-1].shape[1] !=0: # get info gene after SNP (if any); second conditional to evade empty chr
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])] # -1 necessary bcs list of df 0-based; gene_id 0-based by we want following
                mut_annotations['gene2'] = p_anno.loc['gene']
                mut_annotations['locustag2'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag2'] = p_anno.loc['orthologtag']
                mut_annotations['product2'] = p_anno.loc['product']
                mut_annotations['distance2'] = p_anno.loc['loc1'] - mut_annotations['pos'] +1 # +1 to get correct bp distance
                if p_anno.loc['strand'] == 1:
                    mut_annotations['distance2'] = mut_annotations['distance2'] * -1
            # append_annotations intragenic
            #print(mut_annotations['distance1'])
            if ( 'distance1' in mut_annotations and mut_annotations['distance1'] > (-1*promotersize) and mut_annotations['distance1'] < 0) or ( 'distance2' in mut_annotations and mut_annotations['distance2'] > (-1*promotersize) and mut_annotations['distance2'] < 0):
                mut_annotations['type'] = 'P'
            else:
                mut_annotations['type'] = 'I'
            # define nts (repeat of intragenic code)
            if custom_nt_calls_goodpos[i] != 4: # only if ref defined; [0] ok...see next line comment
                mut_annotations['nt_custom'] = nts[custom_nt_calls_goodpos[i]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                mut_annotations['nts'] = mut_annotations['nt_custom']
            # extract derived genotype(s) across all samples
            if hasmutation_goodpos[i] == True:
                if calls_single_sample_goodpos[i] < 4:
                    mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_single_sample_goodpos[i]]
                    mut_annotations['nt_alt'] = nts[calls_single_sample_goodpos[i]]
            # remove duplicates
            mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
        lod_mutAnno.append(mut_annotations)
    dataframe_mut = pd.DataFrame(lod_mutAnno)
    return dataframe_mut

# helper functions
def genomic_position_all(anno_genes_ls,genomelength,chrstarts):
    # get initial placeholder for all genes (differentiate by chr and genic/nongenic)  
    # UPDATE: gene_num > 0 means CDS. changed 'tagnumber' assignment in parse_gff()
    locus_tag_nums = np.ones(genomelength,dtype=float)*0.5 # genome-wide unique CDS tag ('tagnumber'), and intergenic 0.5
    cds_indices = np.ones(genomelength,dtype=float)*0.5 # per chromosome unique; intergenic #chr.5-1; intragenic #gene
    
    for i,this_chr_df in enumerate(anno_genes_ls):
        if this_chr_df.empty: #Skip contigs w/o CDS annotation.
            continue
        
        gene_num = this_chr_df[['tagnumber']].values.flatten() # tdl called cds_num but actually it was not always cds. e.g tRNA; .values.flatten() #1 turn pd.DF to 2D numpy array; #2 turn 2d to 1D 
        genestarts = chrstarts[i] + this_chr_df[['loc1']].values.flatten() ; # .values.flatten() #1 turn pd.DF to 2D numpy array; #2 turn 2d to 1D
        geneends = chrstarts[i] + this_chr_df[['loc2']].values.flatten() ; # .values.flatten() #1 turn pd.DF to 2D numpy array; #2 turn 2d to 1D
        
        for j in range(len(genestarts)-1): # loop over all but last element
            locus_tag_nums[ (genestarts[j]-1):geneends[j] ] = gene_num[j]; # set gene number. Consecutive numbering across all chromosomes; intergenic regions left 0.5

            cds_indices[ (genestarts[j]-1):geneends[j] ] = j+1; # j zero based,but represents chr counter (1-based)
            cds_indices[ (geneends[j]):(genestarts[j+1]-1) ] = j+1+0.5 # intergenic are #chr+0.5; starts with 0.5 prior chr1 which is set when array created
        # mark positions for last gene of chr
        locus_tag_nums[ (genestarts[-1]-1):geneends[-1] ] = gene_num[-1];
        cds_indices[ (genestarts[-1]-1):geneends[-1] ] = len(gene_num); # j zero based,but represents chr counter (1-based)
        # mark trailing cds_intergenic from last gene till end of chromosome (if last chr >> till end of genome)
        if ((i+1) < len(chrstarts)):
            cds_indices[ geneends[-1]:chrstarts[i+1] ] = len(gene_num)+0.5;
        else:
            cds_indices[ geneends[-1]:genomelength ] = len(gene_num)+0.5;
    return [locus_tag_nums,cds_indices]

###################################################################################################
## # # # # # # # # # # # # # All Parallel Evolution Functions # # # # # # # # # # # # # # # # # # #
###################################################################################################

def parallel_evo_module_production_pestis(goodpos2use,contig_positions,annotation_mutations, annotation_genes, params_dict ,plot=True,ortholog=True,plot_title=None, count_multiply_mut_pos=True, dn_ds_only_cand_genes=True,return_prob_sim_cand_genes=False, testing=False, quiet=False,add_muts=False, record_old_locustag=False):
    '''
    Changelog:
    06/2022 IL: added option to not parse orthologs and still output probabilities and plots
    06/2022 IL: added optional plotting output
    09/2022 IL: added optional parsing of 'output_name' to params dictionary to name plot output
    09/2022 IL: added optional plot_title parameter for having plot title
    09/2022 IL: added dynamic xtick plotting
    10/2022 IL: added tight layout for plot export
    11/2022 IL: added parsing of multiple mutation entry in annotation_mutations input (num_mutational_events)
    11/2022 IL: moved dN/dS calculation to helper function, added option to calc dN/dS on all genes or just cand genes (default)
    04/2023 IL: retained only poisson distribution pvalue calculation for significance of num mutations in genes
    06/2024 IL: production outputting of plots, tables

    ## Module to calculate parallel evolution 
    ### Test excess of genes with multiple mutations (thresholds for candidates defined in parameters dict ), recommendation is to use 2 to find all parallely mutated genes above threshold
    ### Test excess of NonSyn (dNdS) in candidate genes compared to expectation given reference genome 
    ### NOTE: dNdS uses simple substitution model (equal P per mut) unless given in params_dict, which should be updated based on observations '''
    #HACKATHON -- decide what to do
    if not quiet: print('Parallel evolution inference.')
    # Find mutations that are adjacent to each other.
    # True for any mutation that is followed by an adjacent mutation, False if not
    # keep trailing bases of adjacent sets (2+)
    chr_pos_gp = contig_positions[goodpos2use,]
    bool_adjacent_mut = np.full( chr_pos_gp.shape[0] , False, dtype=bool) #mutated_genes_pos = []; # keep track of the positions of each event in the table
    for i in range(chr_pos_gp.shape[0]):
        chr = chr_pos_gp[i,0]
        pos = chr_pos_gp[i,1]
        if (i+1) <= (chr_pos_gp.shape[0]-1): # (chr_pos_gp.shape[0]-1) bcs shape[0] not 0-based
            if chr == chr_pos_gp[i+1,0] and (pos+1 == chr_pos_gp[i+1,1]):
                bool_adjacent_mut[i] = True
        
    annotation_mutations['chr_locustag'] = annotation_mutations['chr'].apply(str)+'_'+annotation_mutations['locustag']
    
    # get info for candidate genes
    mutated_genes = np.zeros(0,dtype=int); # keep track of gene number
    mutated_genes_tally = np.zeros(0,dtype=int); # keep track of how many SNPs there are on this gene
    mutated_genes_lengths = np.zeros(0,dtype=int); # keep track of the length of this gene
    chr_locustags_all = np.zeros(0,dtype=object) # record the chr_locustag
    protein_id_all = np.zeros(0,dtype=object) # record the locustag
    orthologtags_all = np.zeros(0,dtype=object) # record the orthologtags
    gene_tally_dict = {} # key=locustag; value=list(first entry tally, second is len, third is nonsyn, 4th is syn)
    total_genic_mutations = 0
    for i, row in annotation_mutations.iterrows():
        if bool_adjacent_mut[i] == False: # ignore leading adjacent mutations
            gene_num_global = row['gene_num_global']
            chr_locustag = row['chr_locustag']
            protein_id = row['protein_id']
            if gene_num_global != 0.5 and gene_num_global != 0: # non-genic 0.5
                total_genic_mutations+=1
                if 'num_mutational_events' in annotation_mutations and count_multiply_mut_pos: ## mutational events counted
                    if chr_locustag in gene_tally_dict:
                        gene_tally_dict[chr_locustag][0] += row['num_mutational_events']
                        if row['NonSyn']: gene_tally_dict[chr_locustag][2] += row['num_mutational_events']
                        else: gene_tally_dict[chr_locustag][3] += row['num_mutational_events']
                    else:
                        mutated_genes_lengths =  (row['loc2']-row['loc1']+1) # +1 bcs loc1/2 1-based, thus underestimate L by 1bp
                        if row['NonSyn']: 
                            gene_tally_dict[chr_locustag] = [row['num_mutational_events'],mutated_genes_lengths, row['num_mutational_events'],0]
                        else: 
                            gene_tally_dict[chr_locustag] = [row['num_mutational_events'],mutated_genes_lengths, 0, row['num_mutational_events']]                        
                        chr_locustags_all = np.append(chr_locustags_all , chr_locustag)
                        protein_id_all = np.append(protein_id_all,protein_id)
                        if ortholog:
                            orthologtags_all = np.append(orthologtags_all , row['orthologtag'])
                        if testing:
                            print(mutated_genes_tally,mutated_genes)
                else:
                    if chr_locustag in gene_tally_dict:
                        gene_tally_dict[chr_locustag][0] += 1 
                        if row['NonSyn']: gene_tally_dict[chr_locustag][2] += 1
                        else: gene_tally_dict[chr_locustag][3] += 1
                    else:
                        mutated_genes_lengths =  (row['loc2']-row['loc1']+1) # +1 bcs loc1/2 1-based, thus underestimate L by 1bp
                        if row['NonSyn']: 
                            gene_tally_dict[chr_locustag] = [1,mutated_genes_lengths, 1 ,0]
                        else: 
                            gene_tally_dict[chr_locustag] = [1,mutated_genes_lengths, 0, 1]                        
                        chr_locustags_all = np.append(chr_locustags_all , chr_locustag)
                        protein_id_all = np.append(protein_id_all,protein_id)
                        if ortholog:
                            orthologtags_all = np.append(orthologtags_all , row['orthologtag'])
                        if testing:
                            print(mutated_genes_tally,mutated_genes)
    mutated_genes_tally=np.array([tally_len[0] for tally_len in gene_tally_dict.values()])
    mutated_genes_lengths=np.array([tally_len[1] for tally_len in gene_tally_dict.values()])
    mutated_genes_nonsyn=np.array([tally_len[2] for tally_len in gene_tally_dict.values()])
    mutated_genes_syn=np.array([tally_len[3] for tally_len in gene_tally_dict.values()])
    mutated_genes=np.array([gene_locustag for gene_locustag in gene_tally_dict.keys()])
    mutated_genes_tally_perGeneLen = mutated_genes_tally/mutated_genes_lengths;
    ## get all gene lengths
    total_gene_lengths=0
    for genomic_element in annotation_genes:
        for i,gene in genomic_element.iterrows():
            total_gene_lengths+=gene['loc2']-gene['loc1']

    #% Define candidates for selection
    mutation_number_threshold = params_dict['Min_num_mutations_cand']; # minimum number of mutations per gene
    mutation_density_threshold = params_dict['Min_mutation_density_cand']; # minimum number of mutations per 1000 bp
    
    mutated_genes_of_interest = ( mutated_genes_tally >= mutation_number_threshold) & (mutated_genes_tally_perGeneLen >= mutation_density_threshold );
    num_mutated_genes_of_interest = sum(mutated_genes_of_interest);
    if not quiet:
        print('Number of genes with multiple mutations: ' + str(num_mutated_genes_of_interest))
        print('Total genic mutations: '+ str(total_genic_mutations))
    output_dir=params_dict['analsysis_params_output_name_folder']
    subprocess.run([f"mkdir -p parallel_evo_results/{output_dir}"], shell=True)
    results_path=f"parallel_evo_results/{output_dir}/parallel_evo_module_results_" + params_dict['output_name'] + ".txt"
    with open(results_path, "w") as myfile:
        myfile.write('#Parallel Evolution output' + "\n")
        myfile.write('#' + params_dict['output_name'] + "\n")
        myfile.write('#Total genic mutations: '+ str(total_genic_mutations) + "\n")
        myfile.write('#Number of genes with multiple mutations: ' + str(num_mutated_genes_of_interest) + "\n")
        myfile.write('#Minimum number mutations required: ' + str(mutation_number_threshold) + "\n")
        myfile.write('#Minimum SNP density for candidate genes: ' + str(mutation_density_threshold) + "\n")        
    
    # break if no candidate genes present
    if num_mutated_genes_of_interest == 0:
        if not dn_ds_only_cand_genes:
            if not quiet:
                print('Calculating dN/dS on all genes, mutations')
                calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations)
            else:
                calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations, print_output=False)
        if not quiet:
            print('No genes with multiple mutation found! >> skip adaptive evolution analysis')
        return [np.array([]),pd.DataFrame(),np.NaN] # return empty array and dataframe
    # get annotation_mutation for SNPs with signature of parallel evolution
    annotation_mutation_paraSignal = annotation_mutations.loc[ annotation_mutations['chr_locustag'].isin( mutated_genes[mutated_genes_of_interest] ) ]

    # =============================================================================
    # Test candidate genes for parallel mutation enrichment
    # =============================================================================
    # NOTE: Simulation based on random genic SNPs in genome and does not use a specific substitution model
    number_of_mutations_on_genome = len(goodpos2use)
    number_of_mutations_genic = int(sum(mutated_genes_tally)); # only genic/CDS!

    [expectedNumberGenesMultipleMutations, expectedNumberOfGenesWithNmutations] = parallel_evolution_counting_and_simulation(number_of_mutations_on_genome,number_of_mutations_genic,params_dict['Min_num_mutations_cand'],params_dict['Min_mutation_density_cand'],params_dict['NumTrialsSim'],params_dict['max_muts_per_gene_to_track'],chr_pos_gp , params_dict['ref_genome_folder'],annotation_genes)
    # expectedNumberGenesMultipleMutations: num candidate genes per sim for parallel evolution based on params_dict['Min_num_mutations_cand'],params_dict['Min_mutation_density_cand']
    # expectedNumberOfGenesWithNmutations: row==sim; col==genes with mutations. col[0] == 1 mutation!, col[1]==2mut...
    
    simProbForObsCand = 1-sum(expectedNumberGenesMultipleMutations < num_mutated_genes_of_interest)/len(expectedNumberGenesMultipleMutations)
    if not quiet:
        print('Simulation-based probability to observe ' + str(num_mutated_genes_of_interest) + ' candidate genes for parallel evolution is: '+str(simProbForObsCand))
        print('Simulated a mean number of candidate genes is ' + str(np.mean(expectedNumberGenesMultipleMutations)))
    with open(results_path, "a") as myfile:
        myfile.write('#Simulation-based probability to observe ' + str(num_mutated_genes_of_interest) + ' candidate genes for parallel evolution is: '+str(simProbForObsCand) + "\n")
        myfile.write('#Simulated a mean number of candidate genes is ' + str(np.mean(expectedNumberGenesMultipleMutations)) + "\n")

    # calc prob to observe candidate mut counts
    locustags_cand = chr_locustags_all[mutated_genes_of_interest]
    protein_id_cand = protein_id_all[mutated_genes_of_interest]
    gene_lengths_cand = mutated_genes_lengths[mutated_genes_of_interest].astype(int)
    if ortholog: orthologtags_cand = orthologtags_all[mutated_genes_of_interest]
    mut_cand_tally = mutated_genes_tally[mutated_genes_of_interest]
    prob_cand_nummut_poisson = np.ones(mut_cand_tally.shape)
    syn_cand=mutated_genes_syn[mutated_genes_of_interest].astype(int)
    nonsyn_cand=mutated_genes_nonsyn[mutated_genes_of_interest].astype(int)
    aa_seq_cand=np.zeros(0,dtype=object) # record the aa sequence
    nt_seq_cand=np.zeros(0,dtype=object) # record the nt sequence
    gene_names_cand=np.zeros(0,dtype=object) # record the gene name
    #genloc_cand=np.zeros(0,dtype=object) # record the genomic_location
    #gene_name_old_locus_tag_cand=np.zeros(0,dtype=object) # record the entry of old locus tag name
    for i,nummut in enumerate(mut_cand_tally):
        this_locus=locustags_cand[i]
        #poisson
        p=number_of_mutations_genic/total_gene_lengths
        prob_cand_nummut_poisson[i] = stats.poisson.sf(nummut,mutated_genes_lengths[mutated_genes_of_interest][i]*p)
        prob_cand_nummut_poisson=prob_cand_nummut_poisson.astype('float64')
        
        # get gene name (when available) 
        gene_names_cand = np.append(gene_names_cand,annotation_mutation_paraSignal.query(f'chr_locustag == "{this_locus}"').iloc[0].gene)
        # produce chrom:start-end (genomic location) (NOTE: would require updating parse_gff() which I am not doing right now)
        # genloc_cand
        # get AA seq from annotation mutations
        aa_seq_cand = np.append(aa_seq_cand,str(annotation_mutation_paraSignal.query(f'chr_locustag == "{this_locus}"').iloc[0].translation))
        # get nt seq from annotation mutations
        nt_seq_cand = np.append(nt_seq_cand,str(annotation_mutation_paraSignal.query(f'chr_locustag == "{this_locus}"').iloc[0].sequence))

    _,prob_cand_nummut_fdr_corrected,_,_=multipletests(prob_cand_nummut_poisson,0.05,'fdr_bh')

    # =============================================================================
    # Generating plot
    # =============================================================================
    if plot:
        plt.rcParams.update({'font.size': 14}) # all label size
        f = figure()
        ax=f.add_subplot(111)
    
        # histogram depicting simulated distribution of expected candidate genes, filter criteria for candidates, P for observation cand count
        mybin = np.arange(max(np.unique(expectedNumberGenesMultipleMutations))+2)-0.5 # +2 in order to now loose last bin due to arange; -0.5 needed for bars align at center with x axis ticks
        plt.hist(expectedNumberGenesMultipleMutations,bins=mybin,rwidth=0.8,color='black',alpha=0.8,density=True)
        max_on_x_axis=int(np.max(np.concatenate([np.array([len(mut_cand_tally)]),np.array(expectedNumberGenesMultipleMutations)])))
        if max_on_x_axis>800:
            step=100
        elif max_on_x_axis>500:
            step=50
        elif max_on_x_axis>100:
            step=20
        elif max_on_x_axis>50:
            step=10
        elif max_on_x_axis>20:
            step=5
        else:
            step=1
        xlim_high=int(max([len(mut_cand_tally),max(np.unique(expectedNumberGenesMultipleMutations))])+1)
        plt.xlim(0,xlim_high)
        plt.xticks(np.arange(0,int(max([len(mut_cand_tally),max(np.unique(expectedNumberGenesMultipleMutations))])+1)+step,step))
        #plt.ylabel('Simulated counts, N='+str(params_dict['NumTrialsSim']))
        plt.ylabel('Proportion of simulations')  
        plt.xlabel('Number of genes with ' + str(params_dict['Min_num_mutations_cand']) + ' or more mutations')
        if plot_title:
            plt.title(plot_title)
        plt.axvline(x=len(mut_cand_tally),color='#e31a1c') 
        # figure out where to plot "observation" text next to axvline
        x_axis_text_plotting_proportion=len(mut_cand_tally)/xlim_high
        if x_axis_text_plotting_proportion>0.5:
            alignment = 'right'
            x_axis_placement=len(mut_cand_tally)-0.025*len(mut_cand_tally)
        else:
            alignment = 'left'
            x_axis_placement=len(mut_cand_tally)+0.025*len(mut_cand_tally)
        updated_ylim=ax.get_ylim()[1]*1.05
        ax.set_ylim(0,updated_ylim)
        text(x_axis_placement,updated_ylim*0.95, "Observation",color='#e31a1c', fontsize=12,horizontalalignment=alignment,verticalalignment='center')
        #text(0.98, 0.95, "min #mut:"+str(params_dict['Min_num_mutations_cand'])+"; min density:"+str(params_dict['Min_mutation_density_cand']), fontsize=12,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
        #text(0.98, 0.88, "P("+ str(len(mut_cand_tally)) + ") = " + str(np.around(simProbForObsCand,3)), fontsize=12,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
        #save pdf
        if 'analsysis_params_output_name_folder' in params_dict:
            output_folder=params_dict['analsysis_params_output_name_folder']+'/'
        else: output_folder=''
        if 'output_name' in params_dict:
            subprocess.run([f"mkdir -p pdf/adaptive_evo/{output_folder} "],shell=True)
            f.savefig(f'pdf/adaptive_evo/{output_folder}' + params_dict['output_name'] + "_" + params_dict['subjectID'] + ".pdf",bbox_inches="tight")
            plt.close()
            if not quiet:
                print(f'Plotted: pdf/adaptive_evo/{output_folder}' + params_dict['output_name'] + "_" + params_dict['subjectID'] + ".pdf")
            with open(results_path, "a") as myfile:
                myfile.write(f'#Plotted: pdf/adaptive_evo/{output_folder}' + params_dict['output_name'] + "_" + params_dict['subjectID'] + ".pdf\n")
        else:
            subprocess.run([f"mkdir -p pdf/adaptive_evo/{output_folder} "],shell=True)
            f.savefig(f'pdf/adaptive_evo/{output_folder}' + params_dict['timestamp'] + "_" + params_dict['subjectID'] + ".pdf",bbox_inches="tight")
            plt.close()
            if not quiet: print(f'Plotted: pdf/adaptive_evo/{output_folder}' + params_dict['timestamp'] + "_" + params_dict['subjectID'] + ".pdf")
            with open(results_path, "a") as myfile:
                myfile.write(f'#Plotted: pdf/adaptive_evo/{output_folder}' + params_dict['timestamp'] + "_" + params_dict['subjectID'] + ".pdf\n")

    # =============================================================================
    # dN/dS calculation + outputting
    # =============================================================================
    # calc if we observe more NS than expected (simulated),
    if dn_ds_only_cand_genes:
        if not quiet:
            print('Calculating dN/dS on only candidate genes')
            calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations, mutated_genes, mutated_genes_of_interest,output_name=results_path)
        else:
            calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations, mutated_genes, mutated_genes_of_interest,print_output=False,output_name=results_path)
    else:
        if not quiet:
            print('Calculating dN/dS on all genes, mutations')
            calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations,output_name=results_path)
        else:
            calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations,print_output=False,output_name=results_path)
    # =============================================================================
    # Generating results tables + outputting
    # =============================================================================
    # calc if we observe more NS than expected (simulated),
    if ortholog:
        res_cand_nummut = np.vstack((orthologtags_cand,locustags_cand,mut_cand_tally,nonsyn_cand,syn_cand,gene_lengths_cand, prob_cand_nummut_poisson)).T
        if not quiet:
            print('Orthologtag Locustag NumMuts NumMutsNonsyn NumMutsSyn GeneLen ProbPoisson')
            print(res_cand_nummut)    
        with open(results_path, "a") as myfile:
            myfile.write('Orthologtag Locustag NumMuts NumMutsNonsyn NumMutsSyn GeneLen ProbPoisson'  + "\n")
            np.savetxt(myfile,  res_cand_nummut,fmt="%s")
    else:
        #if add_muts and record_old_locustag:
        #    res_cand_nummut = np.vstack((protein_id_cand, gene_names_cand,old_locustag_cand,mut_cand_tally,mutations_cand,nonsyn_cand,syn_cand,gene_lengths_cand,prob_cand_nummut_poisson,prob_cand_nummut_fdr_corrected, aa_seq_cand, nt_seq_cand)).T
        #    header='GenomLoc GeneID GeneName GeneName_Locustag NumMuts Mutations NumMutsNonsyn NumMutsSyn GeneLen ProbPoisson FDR_ProbPoisson AASeq NucSeq'.replace(' ','\t')
        #elif record_old_locustag:
        #    res_cand_nummut = np.vstack((protein_id_cand, gene_names_cand,old_locustag_cand,mut_cand_tally,nonsyn_cand,syn_cand,gene_lengths_cand,prob_cand_nummut_poisson,prob_cand_nummut_fdr_corrected, aa_seq_cand, nt_seq_cand)).T
        #    header='GenomLoc GeneID GeneName GeneName_Locustag NumMuts NumMutsNonsyn NumMutsSyn GeneLen ProbPoisson FDR_ProbPoisson AASeq NucSeq'.replace(' ','\t')
        #elif add_muts:
        #    res_cand_nummut = np.vstack((protein_id_cand, gene_names_cand,mut_cand_tally,mutations_cand,nonsyn_cand,syn_cand,gene_lengths_cand,prob_cand_nummut_poisson,prob_cand_nummut_fdr_corrected, aa_seq_cand, nt_seq_cand)).T
        #    header='GenomLoc GeneID GeneName NumMuts Mutations NumMutsNonsyn NumMutsSyn GeneLen ProbPoisson FDR_ProbPoisson AASeq NucSeq'.replace(' ','\t')
        #else: 
        #    res_cand_nummut = np.vstack((protein_id_cand, gene_names_cand,mut_cand_tally,nonsyn_cand,syn_cand,gene_lengths_cand,prob_cand_nummut_poisson,prob_cand_nummut_fdr_corrected, aa_seq_cand, nt_seq_cand)).T
        #    header='GenomLoc GeneID GeneName NumMuts NumMutsNonsyn NumMutsSyn GeneLen ProbPoisson FDR_ProbPoisson AASeq NucSeq'.replace(' ','\t')
        res_cand_nummut = np.vstack((protein_id_cand, gene_names_cand,mut_cand_tally,nonsyn_cand,syn_cand,gene_lengths_cand,prob_cand_nummut_poisson,prob_cand_nummut_fdr_corrected, aa_seq_cand, nt_seq_cand)).T
        header='GeneID GeneName NumMuts NumMutsNonsyn NumMutsSyn GeneLen ProbPoisson FDR_ProbPoisson AASeq NucSeq'.replace(' ','\t')
        if not quiet:
            print(header)
            print(res_cand_nummut)
        with open(results_path, "a") as myfile:
            myfile.write(header + "\n")
            np.savetxt(myfile,  res_cand_nummut,fmt="%s",delimiter='\t')
    if return_prob_sim_cand_genes:
        return [res_cand_nummut,annotation_mutation_paraSignal,simProbForObsCand]
    return [res_cand_nummut,annotation_mutation_paraSignal]

# helper functions
def parallel_evolution_counting_and_simulation(num_mutations_genome, num_mutations_genic , mutation_number_threshold , mutation_density_threshold , numtrials, max_muts_per_gene_to_track, chr_pos_gp, ref_genome_folder, annotation_genes):
    ## Parallel evolution counting and simulation
    # Inputs: genome, number of expected mutations in the genome
    # Output: number of genes with a mutation density (mutations per length of
    # gene) over some threshold

    # num_mutations_genome=number_of_mutations_on_genome
    # num_mutations_genic=number_of_mutations_genic
    # mutation_number_threshold=params_dict['Min_num_mutations_cand']
    # mutation_density_threshold=params_dict['Min_mutation_density_cand']
    # numtrials=params_dict['NumTrialsSim']
    # chr_pos_gp=chr_pos_gp
    # ref_genome_folder=params_dict['ref_genome_folder']
    


    # Get data structures for analysis
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    [locus_tag,chr_tag] = genomic_position_all(annotation_genes, genomelength, chrstarts) #loc_tag numbers all genes(CDS) across genome (intergenic/rna:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    #    locus_tag >> cds_number_by_position_tags
    #    chr_tag >> cds_number_by_position_indices
    # extract vector genelengths that follows order
    
    annotation_genes_sglDF = pd.concat(annotation_genes,sort=False)
    start_pos_gene = annotation_genes_sglDF['loc1'].values
    start_pos_gene = start_pos_gene[~np.isnan(start_pos_gene)] # remove nan, caused by empty chr
    end_pos_gene = annotation_genes_sglDF['loc2'].values
    end_pos_gene = end_pos_gene[~np.isnan(end_pos_gene)] # remove nan, caused by empty chr
    # gene lengths. incl rna genes (however when used those indexes not present in locus_tag (0.5))
    genelengths = (end_pos_gene - start_pos_gene)
    
    expectedNumberGenesMultipleMutations = np.zeros(numtrials) #initialize vector to store simulation results
    expectedNumberOfGenesWithNmutations = np.zeros((numtrials,max_muts_per_gene_to_track));
    # start sim
    for i in range(numtrials):
        # Pick random positions on the genome to mutate
        # Initially get 10x the number of positions you need --> later filter only for num_mutations_genic mutations
        randpos = randint(1, genomelength, 10*num_mutations_genome , dtype=int)
        # Does NOT assume any mutational spectrum!!!
        
        #TDL has scripts for doing the same thing for operon and pathway
        #levels, but the inputs they are generated on may not be relevant for
        #your genome
    
        # Find out in which genes these mutations occurred
        genenums = locus_tag[randpos];
        genenums_intragenic = genenums[genenums>0.5]; #remove intergenic/rna-gene mutations
        genenums_intragenic = genenums_intragenic[0:num_mutations_genic]; # only take as many as you need (generated extra because some removed due to not being on a gene)  
        # NOTE: may be overestimating the number of gene mutations in genic genes vs total num mutations.....
        expectedNumberIntragenicMutations=len(genenums_intragenic)
        # Get a list of the genes mutated, along with the number of times each
        [sim_unique_genes , sim_mutspergene] = np.unique(genenums_intragenic,return_counts=True)
        
        # This calculates the number of mutations per gene length for each gene
        # on which there were simulated mutations
        sim_mutspergenelength = sim_mutspergene/genelengths[sim_unique_genes.astype(int)-1]; #-1 bcs sim_unique_genes 1-based tagnumber > 1st gene is position 0 in genelengths
        sim_mutspergene_above_threshold=sim_mutspergene[(sim_mutspergenelength > mutation_density_threshold)]

        # The final step finds the number of genes that were mutated multiple 
        # times and above the threshold mutation density (per unit length). 
        # Result saved, indexed by trial number.
        expectedNumberGenesMultipleMutations[i] = sum( (sim_mutspergenelength >= mutation_density_threshold) & (sim_mutspergene >= mutation_number_threshold)); # Arolyn, 2018.12.14 > to >= to match input from function
    
        # The second piece of information this script returns is the number of
        # genes with >m mutations
        for j in range(max_muts_per_gene_to_track):
            expectedNumberOfGenesWithNmutations[i,j] = sum( sim_mutspergene_above_threshold >= j+1 ) # +1 due to 0-based
    return [expectedNumberGenesMultipleMutations, expectedNumberOfGenesWithNmutations ]

def codons_in_genome(annotation_genes,allcodons,genewise=False):
    ''' Get probability for each codon across all CDS annotated in the reference genome (annotation_genes), optional output as a per-gene dictionary of codon proportions (NOTE: corrected later if investigating subset of genes)'''
    # possibility to add a flag in order to restrict analysis to genomic region
    annotation_genes_sglDF = pd.concat(annotation_genes,sort=False)
    annotation_genes_sglDF_CDS = annotation_genes_sglDF.loc[ (annotation_genes_sglDF['type'] ==  'CDS') | (annotation_genes_sglDF['type'] ==  'gene') ]
    # Tally codon occurrences over all proteins in genome
    codonCounts = np.zeros( len(allcodons) ,dtype=int ); # for storing tally of codons
    if not genewise:
        for i, row in annotation_genes_sglDF_CDS.iterrows():
            seq = str(row['sequence'])
            startCodon = row['codon_start']
            startPos = startCodon * 3
            codons_gene = [seq[i:i+3] for i in range(startPos, len(seq), 3)]
            codons_gene = collections.Counter(codons_gene) # builds dictionary with all codons (key) in list and counts (value)
            for i,codon in enumerate(allcodons):
                if codon in codons_gene.keys():
                    codonCounts[i] += codons_gene[codon]
        return codonCounts/sum(codonCounts) # probabilities of codons sorted by allcodons order
    else:
        codon_counts_genewise={}
        for i, row in annotation_genes_sglDF_CDS.iterrows():
            codonCounts = np.zeros( len(allcodons) ,dtype=int )
            locustag = str(row['locustag'])
            seq = str(row['sequence'])
            startCodon = row['codon_start']
            startPos = startCodon * 3
            codons_gene = [seq[i:i+3] for i in range(startPos, len(seq), 3)]
            codons_gene = collections.Counter(codons_gene) # builds dictionary with all codons (key) in list and counts (value)
            for i,codon in enumerate(allcodons):
                if codon in codons_gene.keys():
                    codonCounts[i] += codons_gene[codon]
            codon_counts_genewise[locustag]=codonCounts/sum(codonCounts)
        return codon_counts_genewise # probabilities of codons per gene

def mutation_probability(params_dict,annotation_genes,print_statements=True):
    ''' This script examines the probability of a nonsynonymous mutation on a 
        reference genome given some mutational spectrum. '''
    # based on arolyn's m script: @Feb 2019
    ## Define DNA bases, possible mutations, and codons
    # All possible mutations to be considered:
    allbases = np.array(['A','T','G','C']) #fourDNAbases
    allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG']); # 'AT' denotes A to T
    
    # All standard codons
    standard_codon_table = CodonTable.unambiguous_dna_by_id[1]
    allcodons = np.array([c for c in standard_codon_table.forward_table.keys()],dtype=object) # standard codon table, but miss stop codosn
    allcodons = np.append(allcodons , np.array([c for c in standard_codon_table.stop_codons],dtype=object) ) # 64 codons
    # build a combined dictionary (codon>>AA) containing regular and stop codon AA (stop:*)
    codon_all_dict = {}
    for c in allcodons:
        if c in standard_codon_table.forward_table.keys():
            codon_all_dict[c] = standard_codon_table.forward_table[c]
        else:
            codon_all_dict[c] = "*"

    # Generate table of codon composition by base
    codoncompositiontable = codon_composition_table( allbases, allcodons );
    # Rows (64) = all possible codons
    # Columns (4) = number of A's/T's/G's/C's in each codon
    
    ## Generate table of probabilities of nonsynonymous mutations
    # Rows (64) = all possible codons
    # Columns (12) = all possible mutations
    # Entries = probability of nonsynonymous mutation (given a mutation and a codon)
    # Note: If the given mutation cannot be applied to a given codon, the entry is nan. 

    codonnonsyntable = codon_mutation_table( allmuts, allcodons , codon_all_dict );

    ## Calculate mutational spectrum
    # Assuming all types of mutations are equally likely to occur:
    # Get mutational spectrum from experiments, but ignore if fwd or rev strand
    # allmuts = ['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG']
    #AT, TA  0
    #AC, TG  1
    #AG, TC  2
    #GC, CG  3
    #GT, CA  4
    #GA, CT  5
    if params_dict['substitution_spectrum'] == None:
        mutationalspectrum = [1/12] * 12 # uniform distribution. replace with time stratified observation based on all samples
    else:
        afile = open(params_dict['substitution_spectrum'], 'rb')
        mutationalspectrum = pickle.load(afile)
        afile.close()
        if print_statements: print('Observed substitution spectrum loaded')

    ## Calculate codon distribution in reference genome ordered as in allcodons
    codondistribution = codons_in_genome( annotation_genes , allcodons, False );
    
    ## Calculate probability of nonsynonymous mutation
    # Takes into account: mutation spectrum, abundance of codons on genome,
    # abundance of bases in codons, probability of a given mutation being
    # nonsynonymous on a givne codon...
    probnonsyn = 0; # probability of nonsynonymous mutations over all possible mutations

    for i,mut in enumerate(allmuts): # loop through all possible mutations
        # Probability that this mutation occurs
        prob_base_base = mutationalspectrum[i]
        # Find the codons that can undergo this mutation:
        base = mut[0] # base that gets mutated; ex. A
        baseindex = np.where(allbases==base) # ex. A is indexed in position 1
        basecodonoccurrences = codoncompositiontable[:,baseindex].flatten()  # ex. how many A's in each codon
        # Indices of codons that have the relevant initial base:
        basecodonoccurrences_bool = basecodonoccurrences > 0; # ex. AAT has A's but GGC does not
        # Probability that this mutation occurs on a given relevant codon
        # Take into account base composition of codons
        basecountincodon = basecodonoccurrences[ basecodonoccurrences_bool ];
        # Take into account codon abundance on reference genome
        probcodonongenome = codondistribution[ basecodonoccurrences_bool ]
        # Combine these two probabilities
        probmutoncodon = basecountincodon*probcodonongenome
        # Renormalize (sum = 1 over all relevant codons)
        probmutoncodon = probmutoncodon/sum(probmutoncodon);

        # Probability that this mutation is nonsynonymous at each relevant codon
        thismutnonsynoncodon = codonnonsyntable[:,i];
        probmutnonsynoncodon = thismutnonsynoncodon[basecodonoccurrences_bool];
        # Overall probability that this mutation is nonsynonymous over all possible codons
        probmutnonsyn=prob_base_base*sum(probmutoncodon*probmutnonsynoncodon);

        # Add contribution of this mutation to the total probability of a nonsynonymous mutation
        probnonsyn += probmutnonsyn  
        
    if print_statements: print('Probability of nonsynonymous mutation across genome: ' + str(probnonsyn) )
    return probnonsyn # Prob for N occuring 

def codon_composition_table( allbases, allcodons ):
    ''' Build table containing base counts (col) for each codon (row) '''
    codoncompositiontable = np.zeros((len(allcodons),len(allbases)),dtype=int)
    for i,codon in enumerate(allcodons):
        for j,base in enumerate(allbases):
            codoncompositiontable[i,j] = codon.count(base)
    return codoncompositiontable

def codon_mutation_table(allmuts , allcodons , codon_all_dict):
    ''' Generates table of probabilities that a given mutation on a codon is nonsynonymous '''    
    table = np.zeros( (len(allcodons), len(allmuts) ) ); # Initialize table
    for i,codon in enumerate(allcodons):
        for j,mut in enumerate(allmuts):
            # Calculates the probability that a given mutation is nonsynonymous on a given codon and then updates the table
            table[i,j] = prob_nonsyn_codon_mutation( codon, mut , codon_all_dict);
    return table

def prob_nonsyn_codon_mutation(codon,mut,codon_all_dict):
    ''' Calculates the probability that a given mutation leads to a nonsynonymous change across a given codon '''
    aa0 = codon_all_dict[codon] # AA of codon
            
    # Find the positions on the codon at which mutation could occur
    possiblemuts=[i for (i, base) in enumerate(codon) if base == mut[0]]
    ctr = 0
    if len(possiblemuts) == 0: # if the mutation cannot occur on this codon
        probability = float('nan');
    else: # mut can occur at least once
        for pos in possiblemuts:
            newcodon = list(codon)
            newcodon[pos] = mut[1] # mutate codon position that carries mut[0]
            aa1 = codon_all_dict[ "".join(newcodon) ]
            if aa0 != aa1:
                ctr += 1
        probability = ctr/len(possiblemuts) # fraction of mutations that were nonsynonymous
    return probability
     
def calculate_genome_wise_dn_ds(params_dict,annotation_genes,annotation_mutations_from_para, mutated_genes=[], genes_of_interest=[], save_results=True, print_output=True,return_values=False, calc_proportion_confint=False, quiet = False,output_name=''):
    """
    Function to calculate dN/dS ratio and compare to theoretical expectations at a genome wide level
    Set of total mutations can be subset to specific genes providing a set of chr_locustags for annotation_mutations
    
    was previously implemented in parallel evo module, which now calls this as a helper function
    IL 14.02.2023 added optional outputting of probs and CI, rather than just printing/saving results
    """
    if print_output: print('Investigate dN/dS')
    if output_name == '':
        output_file="dn_ds_results/dn_ds_results_" + params_dict['output_name'] + ".txt"
    else:
        output_file = output_name
    ## Gene numbers of genes of interest, if provided use subset, if not use all mutations (genic filtered later during counting)
    if len(genes_of_interest) > 0 and len(mutated_genes) > 0:
        cand_genes_chr_locustags = mutated_genes[genes_of_interest]; # these_gene_nums
    else:
        cand_genes_chr_locustags = annotation_mutations_from_para['chr_locustag'].unique()

    ## check if num_mutational_events is present
    if 'num_mutational_events' in annotation_mutations_from_para.columns:
        use_multiple_events_count = True
    else: use_multiple_events_count = False
    # Observed dNdS
    cand_muts_N = 0;
    cand_muts_S = 0;
    for chr_locustag in cand_genes_chr_locustags:
        cand_mut_anno = annotation_mutations_from_para.loc[ annotation_mutations_from_para['chr_locustag'] == chr_locustag ]
        for i,row in cand_mut_anno.iterrows():
            ## ensures only count if in genic regions
            if row['NonSyn']:
                if use_multiple_events_count:
                    cand_muts_N += row['num_mutational_events']
                else:
                    cand_muts_N += 1
            else:
                if use_multiple_events_count:
                    cand_muts_S += row['num_mutational_events']
                else:
                    cand_muts_S += 1
    if print_output:
        if (cand_muts_N+cand_muts_S) == 0: # necessary in case no S/N mutations but all P/U
            print("No dN/dS calculated >> no S/N mutations.")
        else:
            probN_obs = cand_muts_N/(cand_muts_N+cand_muts_S)        
            # Simulate genome specific expected probability for NonSyn to occur
            # simulate only if value not yet calculated. > stored in regenome/probNsim.pk     
            # probN_sim = 0.6829; 
            probN_sim = mutation_probability(params_dict,annotation_genes,print_output)
            # NOTE: mutation_probability uses a pre-calculted substitution spectrum (mutationalspectrum) when specified in params or a uniform mutation spectrum.
            
            # Binomial test
            prob_obsN_excess = stats.binom_test(cand_muts_N, n=(cand_muts_N+cand_muts_S), p=probN_sim, alternative='greater')
            prob_obsS_excess = stats.binom_test(cand_muts_S, n=(cand_muts_N+cand_muts_S), p=(1-probN_sim), alternative='greater')
            print('Observed P(N): ' + str(probN_obs) + '; N: ' + str(cand_muts_N) + ', S: ' + str(cand_muts_S))
            print('Expected P(N): ' + str(probN_sim) )
            print('Probability of enrichment in observed count of N: ' + str(prob_obsN_excess) )
            print('Probability of enrichment in observed count of S: ' + str(prob_obsS_excess) )
            if save_results:
                subprocess.run([f"mkdir -p dn_ds_results"], shell=True)
                with open(output_file, "a") as myfile:
                    myfile.write('#Observed P(N): ' + str(probN_obs) + '; N: ' + str(cand_muts_N) + ', S: ' + str(cand_muts_S) + "\n")
                    myfile.write('#Expected P(N): ' + str(probN_sim) + "\n")
                    myfile.write('#Probability of enrichment in observed count of N: ' + str(prob_obsN_excess) + "\n")
                    myfile.write('#Probability of enrichment in observed count of S: ' + str(prob_obsS_excess) + "\n")
                print(f'Saved results in file: {output_file}')
    else:
        if (cand_muts_N+cand_muts_S) == 0:
            return (np.NaN,np.NaN)
        probN_obs = cand_muts_N/(cand_muts_N+cand_muts_S)  
        probN_sim = mutation_probability(params_dict,annotation_genes,print_output)
        # NOTE: mutation_probability uses a pre-calculted substitution spectrum (mutationalspectrum) when specified in params or a uniform mutation spectrum.
        # Binomial test
        prob_obsN_excess = stats.binom_test(cand_muts_N, n=(cand_muts_N+cand_muts_S), p=probN_sim, alternative='greater')
        prob_obsS_excess = stats.binom_test(cand_muts_S, n=(cand_muts_N+cand_muts_S), p=(1-probN_sim), alternative='greater')
        if save_results:
            subprocess.run([f"mkdir -p dn_ds_results"], shell=True)
            with open(output_file, "a") as myfile:
                myfile.write('#Observed P(N): ' + str(probN_obs) + '; N: ' + str(cand_muts_N) + ', S: ' + str(cand_muts_S) + "\n")
                myfile.write('#Expected P(N): ' + str(probN_sim) + "\n")
                myfile.write('#Probability of enrichment in observed count of N: ' + str(prob_obsN_excess) + "\n")
                myfile.write('#Probability of enrichment in observed count of S: ' + str(prob_obsS_excess) + "\n")
        if return_values==False:
            if calc_proportion_confint:
                N_prop_conf=proportion_confint(cand_muts_N,cand_muts_N+cand_muts_S,0.05,method='wilson')
                S_prop_conf=proportion_confint(cand_muts_S,cand_muts_N+cand_muts_S,0.05,method='wilson')
                return (prob_obsN_excess,prob_obsS_excess,N_prop_conf,S_prop_conf)
            return (prob_obsN_excess,prob_obsS_excess)
        else:
            return (prob_obsN_excess,prob_obsS_excess,probN_obs)

def calculate_dn_ds_across_genes(params_dict,annotation_genes,annotation_mutation,calc_proportion_confint=False,precalculated_gene_dict=None,add_one_numer_and_denom=False):
    """
    Calculates dnds ratio for genes given in annotation_mutations, based upon overall likelihood of a mutation in these genes leading to syn or nonsyn mut

    Optional return of 95% CI for estimate
    """
    if precalculated_gene_dict:
        dn_prob_per_gene=precalculated_gene_dict
    else: dn_prob_per_gene=cacluate_expected_dn_ds_per_gene(params_dict,annotation_genes)
    locustags_investigated=[]
    annotation_genes_sglDF = pd.concat(annotation_genes,sort=False)
    cand_muts_N=0
    cand_muts_S=0
    for locustag in np.unique(annotation_mutation['locustag'][annotation_mutation['type'].isin(['N','S'])]):
        locustags_investigated.append(locustag)
        this_locustag_annotation_mutations=annotation_mutation[annotation_mutation['locustag'].isin([locustag])]
        for i,row in this_locustag_annotation_mutations.iterrows():
            ## ensures only count if in genic regions
            if row['NonSyn']:
                if 'num_mutational_events' in this_locustag_annotation_mutations.columns:
                    cand_muts_N += row['num_mutational_events']
                else:
                    cand_muts_N += 1
            elif row['type'] == "S":
                if 'num_mutational_events' in this_locustag_annotation_mutations.columns:
                    cand_muts_S += row['num_mutational_events']
                else:
                    cand_muts_S += 1
    # get proportion weighting each gene should get for contribution to expected nonsyn muts (weight by proportion of total gene lengths)
    seq_lengths=[]
    for locustag in locustags_investigated:
        seq_lengths.append(len(annotation_genes_sglDF[annotation_genes_sglDF['locustag'].isin([locustag])]['sequence']))
    weighting=np.array(seq_lengths)/np.sum(np.array(seq_lengths))
    overall_prob=0
    for index,locustag in enumerate(locustags_investigated):
        overall_prob+=dn_prob_per_gene[locustag]*weighting[index]
    if add_one_numer_and_denom:
        cand_muts_N+=1
        cand_muts_S+=1
    sum_all_genic_mutations=(cand_muts_N+cand_muts_S)
    if sum_all_genic_mutations >0:
        prob_obsN_excess = stats.binom_test(cand_muts_N, n=sum_all_genic_mutations, p=overall_prob, alternative='greater')
        prob_obsS_excess = stats.binom_test(cand_muts_S, n=sum_all_genic_mutations, p=(1-overall_prob), alternative='greater')
        # calculate binom
        if calc_proportion_confint:
            N_prop_conf=proportion_confint(cand_muts_N,sum_all_genic_mutations,0.05,method='wilson')
            S_prop_conf=proportion_confint(cand_muts_S,sum_all_genic_mutations,0.05,method='wilson')
            return (cand_muts_N/sum_all_genic_mutations,cand_muts_S/sum_all_genic_mutations,N_prop_conf,S_prop_conf,overall_prob,sum_all_genic_mutations)
        else: return (prob_obsN_excess,prob_obsS_excess,overall_prob)
    else: 
        if calc_proportion_confint: return (np.NaN,np.NaN,np.NaN,np.NaN,overall_prob)
        else:
            return (np.NaN,np.NaN,overall_prob)

def cacluate_expected_dn_ds_per_gene(params_dict,annotation_genes):
    ''' This script examines the probability of a nonsynonymous mutation on a 
    reference genome given some mutational spectrum. '''
    # HACKATHON
    # based on arolyn's m script: @Feb 2019
    ## Define DNA bases, possible mutations, and codons
    # All possible mutations to be considered:
    allbases = np.array(['A','T','G','C']) #fourDNAbases
    allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG']); # 'AT' denotes A to T
    
    # All standard codons
    standard_codon_table = CodonTable.unambiguous_dna_by_id[1]
    allcodons = np.array([c for c in standard_codon_table.forward_table.keys()],dtype=object) # standard codon table, but missing stop codons
    allcodons = np.append(allcodons , np.array([c for c in standard_codon_table.stop_codons],dtype=object) ) # 64 codons
    # build a combined dictionary (codon>>AA) containing regular and stop codon AA (stop:*)
    codon_all_dict = {}
    for c in allcodons:
        if c in standard_codon_table.forward_table.keys():
            codon_all_dict[c] = standard_codon_table.forward_table[c]
        else:
            codon_all_dict[c] = "*"

    # Generate table of codon composition by base
    codoncompositiontable = codon_composition_table( allbases, allcodons );
    # Rows (64) = all possible codons
    # Columns (4) = number of A's/T's/G's/C's in each codon
    
    ## Generate table of probabilities of nonsynonymous mutations
    # Rows (64) = all possible codons
    # Columns (12) = all possible mutations
    # Entries = probability of nonsynonymous mutation (given a mutation and a codon)
    # Note: If the given mutation cannot be applied to a given codon, the entry is nan. 

    codonnonsyntable = codon_mutation_table( allmuts, allcodons , codon_all_dict );

    ## Calculate mutational spectrum
    # Assuming all types of mutations are equally likely to occur:
    # Get mutational spectrum from experiments, but ignore if fwd or rev strand
    # allmuts = ['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG']
    #AT, TA  0
    #AC, TG  1
    #AG, TC  2
    #GC, CG  3
    #GT, CA  4
    #GA, CT  5
    if params_dict['substitution_spectrum'] == None:
        mutationalspectrum = [1/12] * 12 # uniform distribution. replace with time stratified observation based on all samples
    else:
        afile = open(params_dict['substitution_spectrum'], 'rb')
        mutationalspectrum = pickle.load(afile)
        afile.close()
    ## Calculate codon distribution in reference genome ordered as in allcodons
    codondistribution = codons_in_genome( annotation_genes , allcodons, True );

    # Calculate probability of nonsynonymous mutations per gene and output as dict for gene locustag --> prob    
    locustag_based_nonsyn_prob={}
    for genelocustag in codondistribution:
        probnonsyn = 0; # probability of nonsynonymous mutations over all possible mutations

        for i,mut in enumerate(allmuts): # loop through all possible mutations
            # Probability that this mutation occurs
            prob_base_base = mutationalspectrum[i]
            # Find the codons that can undergo this mutation:
            base = mut[0] # base that gets mutated; ex. A
            baseindex = np.where(allbases==base) # ex. A is indexed in position 1
            basecodonoccurrences = codoncompositiontable[:,baseindex].flatten()  # ex. how many A's in each codon
            # Indices of codons that have the relevant initial base:
            basecodonoccurrences_bool = basecodonoccurrences > 0; # ex. AAT has A's but GGC does not
            # Probability that this mutation occurs on a given relevant codon
            # Take into account base composition of codons
            basecountincodon = basecodonoccurrences[ basecodonoccurrences_bool ];
            # Take into account codon abundance on reference genome
            probcodonongenome = codondistribution[genelocustag][ basecodonoccurrences_bool ]
            # Combine these two probabilities
            probmutoncodon = basecountincodon*probcodonongenome
            # Renormalize (sum = 1 over all relevant codons)
            probmutoncodon = probmutoncodon/sum(probmutoncodon);

            # Probability that this mutation is nonsynonymous at each relevant codon
            thismutnonsynoncodon = codonnonsyntable[:,i];
            probmutnonsynoncodon = thismutnonsynoncodon[basecodonoccurrences_bool];
            # Overall probability that this mutation is nonsynonymous over all possible codons
            probmutnonsyn=prob_base_base*sum(probmutoncodon*probmutnonsynoncodon);

            # Add contribution of this mutation to the total probability of a nonsynonymous mutation
            probnonsyn += probmutnonsyn  
        locustag_based_nonsyn_prob[genelocustag]=probnonsyn
    return locustag_based_nonsyn_prob # Prob for N occuring 

def get_pairwise_genetic_distance(sample1,sample2,calls,sampleNames):
    """Calculates pairwise genetic distance, normalized to the number of comparisons possible"""
    col_index_to_compare=np.nonzero(np.in1d(sampleNames, sample1))[0]
    row_index_to_compare=np.nonzero(np.in1d(sampleNames, sample2))[0]
    calls_same_for_comparissons=np.sum(((calls[:,col_index_to_compare] == calls[:,row_index_to_compare]) & (calls[:,col_index_to_compare] != 4) & (calls[:,row_index_to_compare] != 4)))
    denominator_for_comparissons=np.sum(((calls[:,col_index_to_compare] !=4 ) & (calls[:,row_index_to_compare] != 4)))
    normalized_genetic_distance=(denominator_for_comparissons-calls_same_for_comparissons)/denominator_for_comparissons
    return normalized_genetic_distance


###################################################################################################
## # # # # # # # # # # # # # All phylogenetic-related analyses # # # # # # # # # # # # # # # # # ##
###################################################################################################
# Tree I/O
def parse_tree(tree_path):
    """helper function to parse tree file, with error catching. Returns first tree in a file, if multiple"""
    try:
        parsed_tree=Phylo.read(tree_path,'newick')
    except ValueError:
        parsed_tree=None
        print("Multiple trees in file, trying to parse first (filled) tree in file.")
        trees= Phylo.parse(tree_path,'newick')
        for tree in trees:
            if tree.count_terminals() > 1:
                parsed_tree=tree
                print("Found filled tree, with length", tree.count_terminals())
                break
    return parsed_tree


# Tree structure-related
def reroot_tree(tree_path, outgroup_name_for_lca = "", verbose = False):
    if verbose: print('Note: Uses first strain as LCA, unless outgroup_name_for_lca provided as a string.')
    if outgroup_name_for_lca != "": print("Outgroup name provided, will attempt to use sample", outgroup_name_for_lca, "as outgroup.")
    rc('font', **{'sans-serif':'Arial'})
    if verbose: print("Attempting to reroot tree, based on provided outgroup_name_for_lca",outgroup_name_for_lca)
    input_tree=Phylo.read(tree_path, "newick")
    input_tree.root_with_outgroup(outgroup_name_for_lca)
    tree_path=tree_path.replace("_latest","_latest_rerooted_"+outgroup_name_for_lca)            
    Phylo.write(input_tree, tree_path, "newick")
    if verbose: 
        print()
        print(f"Successfully rerooted tree to filename",tree_path)

def find_clade_terminals(clade_name, tree):
    """Finds a sample name in a tree and returns all sample names in a list at or below this branch in the clade"""
    try:
        parsed_tree=Phylo.read(tree,'newick')
    except ValueError:
        parsed_tree=None
        print("Multiple trees in file, trying to parse first (filled) tree in file.")
        trees= Phylo.parse(tree,'newick')
        for tree in trees:
            if len(tree.get_terminals()) > 1:
                parsed_tree=tree
                print("Found filled tree, with length", len(tree.get_terminals()))
                break
    clades_found=[]
    for clade in parsed_tree.find_elements(clade_name):
        clades_found.append(clade)
    if len(clades_found) == 0:
        print("Error, no terminal branch found with input clade_name:",clade_name)
    else:
        parent=parsed_tree.get_path(clades_found[0])[-2] ## find parent node
        ancestor_clades=parent.get_terminals()
        output=[]
        for ancestor in ancestor_clades:
            output.append(ancestor.name)
        return output


# ancestral reconstruction
def create_ancestral_reconstruction(treename,treesampleNamesLong_ref_outgroup,calls_for_tree_ref_outgroup, outdir="ancestral_reconstruction"):
    """Creates internal node fasta SNP alignments from input tree, using treetime"""
    subprocess.run([f"rm -fr {outdir} ; mkdir -p {outdir}" ],shell=True)
    with open(f'{outdir}/snp_table_treetime.fasta','w') as file:
        for index,sample in enumerate(treesampleNamesLong_ref_outgroup):
            str_alignment=''.join(calls_for_tree_ref_outgroup[:,index])
            for_output=str_alignment.replace('?','N')
            file.write(f'>{sample}')
            file.write('\n')
            file.write(f'{for_output}')
            file.write('\n')
    # run treetime
    subprocess.run([f"treetime ancestral --aln {outdir}/snp_table_treetime.fasta --tree {treename} --outdir {outdir}"],shell=True)

def parse_ancestra_fasta_treetime(fasta, node_label_to_use=None):
    """
    Loads in ancestral sequence fasta output from treetime. Note this output is only the indices that went into the calls for tree (usually treetime)
    """
    if type(fasta) == str:
        ref=SeqIO.parse(fasta,'fasta')
        for record in ref:
            if node_label_to_use:
                if record.id == node_label_to_use:
                    records = nts2idx(np.array([x for x in record.seq]))
    return records


## Generate trees colored by basecall for each SNP functions
def generate_mutation_colored_tree_helper_generate_temp_tree(tree, samples_basecalls_csv, outgroup_name, color_lca, label_lca, indels, mutational_events_for_name=[]):
    """
    Generates new trees for each mutation (position) with coloring based on base call at the position
    iterates over all the mutations in samples_basecalls_csv (csv generated in build_table_for_tree_labling)
    samples_basecalls_csv is csv of chrom,pos,samples.... on first line, following lines are chrom,pos, SNP call for each sample


    Changelog:
    06.2022 IL: added docstring, parameter options for: 1) designating outgroup, 2) counting mutation events for naming, 3) coloring reference based on basecall, 4) labelling reference name with basecall
    09.2022 IL: added parsing option for indels
    03.2023 IL: added mutational_events_for_name param which is num_mutational_events :)
    """
    if len(mutational_events_for_name)>0:
        print("Note: Output name of tree will include number of mutational events for this mutation (position), based on current tree structure")
    f=open(samples_basecalls_csv,'r').readlines()
    lca_index_found=False ## flag to find outgroup name only once
    for i, line in enumerate(f[1:]):
        ## start iterating over each mutation (position), skip header
        l=line.strip().split(',')
        if len(l) < 5:
            continue
        chromosome=l[0]
        pos=str(int(l[1])) # turn 0-based pos into 1-based pos.
        if not lca_index_found:
            if len(outgroup_name) > 1:
                lca_index=f[0].strip().split(',').index(outgroup_name)
            else:
                #use first strain as lca
                lca_index=2
            lca_index_found=True
            print("lca index =",lca_index, "corresponding to sample name", f[0].strip().split(',')[lca_index])
        lca=l[lca_index]
        #count mutations
        newtree_nexus_snp_colored = generate_mutation_colored_tree_helper_add_snp_call(tree, samples_basecalls_csv, chromosome, pos,outgroup_name, color_lca, label_lca,indels)
        pos_for_output=str(int(pos)+1)
        if len(mutational_events_for_name)>0:
            num_mutational_events=int(mutational_events_for_name[i])
        #if a==0:
            #print('NO MUTS:')
            #print(chromosome, pos)
        #save trees
            f1=open(f'{chromosome}_{pos_for_output}_{num_mutational_events}.tree','w')
        else:
            f1=open(chromosome+'_'+pos_for_output+'.tree','w')
        t=open('tempnexus.txt','r').readlines()
        for line in t:
            f1.write(line)

def generate_mutation_colored_tree_helper_add_snp_call(tree, bascalls_per_sample_csv, chromosome, pos, reference_name_in_tree, color_lca, label_lca,indels):
    """
    creates temptree.txt tsv file w/ basecall for each sample, then calls generate_mutation_colored_tree_helper_add_colors_to_basecall_tree to generate new tree file with basecall as color labels

    changelog:
    6.2022 IL: added docstring, removed limit to sample numbers and naming convention required to add sample names to temptree.txt, renamed dict --> bascalls_per_sample_csv for clarity, added parsing for generate_mutation_colored_tree_helper_add_colors_to_basecall_tree
    """
    f=open(bascalls_per_sample_csv).readlines()
    fo=open('temptree.txt','w')
    header=f[0].strip().split(',')
    locations=[]
    for n,i in enumerate(header):
    #        if (i.startswith('S') or i.startswith('D') or i.startswith('R') or i.startswith('P')) and len(i)<100:
        locations.append(n)
    for line in f:
        l=line.strip('\n').split(',')
        if l[0]==chromosome and l[1]==pos:
            for i in locations:
                if len(l) > i and len(l[i])>0:
                    fo.write(header[i]+'\t'+l[i]+'\n')
                else:
                    fo.write(header[i]+'\t?\n')
                #if i > len(l):
                    #print(line, l, i)
            break
    fo.close()
    newtree_nexus_snp_colored = generate_mutation_colored_tree_helper_add_colors_to_basecall_tree(tree,'temptree.txt',reference_name_in_tree, color_lca, label_lca,indels,chromosome,pos)
    return newtree_nexus_snp_colored

def generate_mutation_colored_tree_helper_add_colors_to_basecall_tree(tree, dictionary, reference_name_in_tree, color_lca, label_lca, indels,chromosome,pos):
    """
    generates tempnexus.txt with color labels and tiplabels based on basecall
    
    Changelog:
    09.2022 IL: Added parsing of indels and coloring based on classes, coloring classes will repeat if >6 indel classes exist (unlikely)
    """
    f=open(dictionary).readlines()
    numStrains=len(f)
    annotation={}
    intree_labels_to_outtree_labels_with_snp={}
    #print header for tempnexus
    fo=open('tempnexus.txt','w')
    fo.write('#NEXUS\nbegin taxa;\n\tdimensions ntax='+str(numStrains-2)+';\n\ttaxlabels\n')
    if not indels:
        colors={'A':'[&!color=#-16776961]', 'C':'[&!color=#-16725916]','G':'[&!color=#-3670016]','T':'[&!color=#-3618816]','?':'[&!color=#-16777216]'}
    else:
        colors={'0':'[&!color=#-10223416]','?':'[&!color=#-16777216]'}
        ## not all colors will be used, repeats to ensure it has enough anyways
        colors_multiple_indels=['[&!color=#-3670016]','[&!color=#-32768]','[&!color=#-3618816]','[&!color=#-65408]','[&!color=#-16725916]','[&!color=#-16776961]']
    ambiguous=['R','Y','M','K','S','W']
    for each in ambiguous:
        colors[each]='[&!color=#-16777216]'
    indels_colored={}
    #get annotations
    f=open(dictionary, 'r').readlines()
    for line in f:
        if not line.startswith('#'):
            l=line.split()
            annotation[l[0]]=l[1]
    # combine names and annotations, dict sample --> sample--*A/T/C/G
    for i in annotation.keys():
        if not label_lca and i == reference_name_in_tree:
            intree_labels_to_outtree_labels_with_snp[i] = i
        else: 
            intree_labels_to_outtree_labels_with_snp[i] = i + '--*' + annotation[i]
    #make new tree
    f=open(tree,'r').readlines()
    for line in f:
        ## first need order of samples as they appear in tree, to put color labels in correct order
        sample_order_in_tree=[]
        for sample in line.strip().replace('(','').replace(')','').split(','): ## each sample will have own string, with other info from previous tree
            sample_name=sample.split(':')[0] ## separate sample_name from info from tree (eg relatedness, not important for us right now)
            sample_order_in_tree.append(intree_labels_to_outtree_labels_with_snp[sample_name]) ## add to list of sample order in tree
            samplename_with_basecall=sample_order_in_tree[-1] ## samplename just added
            if not indels:
                if not color_lca and sample_name == reference_name_in_tree: 
                    fo.write('\t\''+samplename_with_basecall+'\'[&!color=#-16777216]\n')
                else:
                    fo.write('\t\''+samplename_with_basecall+'\''+colors[samplename_with_basecall[-1]]+'\n') #write down new colors to nexus file
            else:
                if not color_lca and sample_name == reference_name_in_tree: 
                    fo.write('\t\''+samplename_with_basecall+'\'[&!color=#-16777216]\n')
                else:
                    if samplename_with_basecall[-1] not in colors:
                        if samplename_with_basecall[-1] in indels_colored:
                            color = indels_colored[samplename_with_basecall[-1]]
                        else:
                            index_to_add=len(indels_colored) % len(colors_multiple_indels)
                            if len(indels_colored) > len(colors_multiple_indels): print(f"Coloring classes of indels will repeat for tree {chromosome}_{pos}")
                            indels_colored[samplename_with_basecall[-1]]=colors_multiple_indels[index_to_add]
                            color = indels_colored[samplename_with_basecall[-1]]
                        fo.write('\t\''+samplename_with_basecall+'\''+color+'\n')
                    else:
                        fo.write('\t\''+samplename_with_basecall+'\''+colors[samplename_with_basecall[-1]]+'\n')
        fo.write(';\nend\n\nbegin trees;\n\ttree tree_1=[&R] ')
        ## next, need to replace all samplenames in tree file with newly updated labels with colors, as above
        for intree_labels,outtree_labels in intree_labels_to_outtree_labels_with_snp.items():
            line=line.replace(intree_labels,outtree_labels)
        newtree_nexus_snp_colored=line
    for line in newtree_nexus_snp_colored:
        fo.write(line)
    fo.write('end;\n\nbegin figtree;\n')
    fo.write('\tset tipLabels.fontSize=10;\n')
    fo.write('end;')
    fo.close()
    return newtree_nexus_snp_colored

def generate_mutation_colored_tree(tree,samples_basecalls_csv, outgroup_name_for_lca="", count_mutations_for_name=[], color_lca=False, label_lca=True, indels=False):
    """
    Main function for generating trees for every SNP
    Calls generate_mutation_colored_tree_helper_generate_temp_tree and sets font. 
    Set values for downstream helper function calls in this function. 
    
    Changelog:
    06.2022 IL: added additional print statements for flagging information, expected output behavior, added additional flags for modifying output
    06.2022 IL: added optional rerooting on all SNP trees to user input LCA, using biopython's Phylo, modified parameter names to lca instead of reference, outputting parameters used for run
    09.2022 IL: added parsing option for indel calls.
    03.2023 IL: externalized rerooting of tree to separate function, to be called prior to coloring tree.
    """
    print("Note: Coloring last common ancestor (lca) tip by basecall option is currently",color_lca, "... set color_lca to True/False to modify")
    print("Note: Labelling last common ancestor with basecall is currently", label_lca, "... set label_lca to True/False to modify")
    # save parameters into output directory
    paramters_df=pd.DataFrame({"parameter_name":list(locals().keys()),"values_input":list(locals().values())})
    paramters_df.to_csv("parameters.tsv",sep="\t",index=False)
    ## print parameters and log for start, end at runtime
    if indels:
        print("Generating per-indel trees, uncalled black, indel as red, no indel call as green.")
    else:
        print("Generating per-snp trees, using basecalls as color.")
    generate_mutation_colored_tree_helper_generate_temp_tree(tree,samples_basecalls_csv,outgroup_name_for_lca, color_lca, label_lca, indels, count_mutations_for_name)
    if indels:
        print("Successfully generated per-indel trees, uncalled black, Reference (no indel) as purple, other indel classes as other colors.")
    else:
        print('Successfully generated per-snp trees, colored by mutation basecall')

def build_table_for_tree_labeling(p_chr_table,treeSampleNamesLong,calls_for_tree,patient="",indels=False):
    """Purges any previous per-SNP trees and recreates folder structure to fill. Creates for_tree_labelling.csv in appropriate locaiton (which is just a SNP table) to color trees in later functions."""
    # make new folder ('tree_counting'), wipe all previous content, add table
    if calls_for_tree.shape[1] != treeSampleNamesLong.shape[0]:
        print("WARNING: inputs for treeSampleNamesLong and calls_for_tree inputs have different dimensions (differing number of samples), downstream errors will occur with indexing")
        print("Aborting function execution.")
        return
    if patient != "":
        subprocess.run(["rm -fr tree_counting/subject"+patient+" ; mkdir tree_counting/subject"+patient ],shell=True)
        with open("tree_counting/subject"+patient+"/for_tree_labeling.csv",'w') as csv_file:
            csv_file.write(",".join(np.append(np.array(['chr','pos']),treeSampleNamesLong))+"\n") # write header
            for i in range(p_chr_table.shape[0]):
                csv_file.write(",".join(np.append( np.array([str(p_chr_table[i,0]),str(p_chr_table[i,1])]) ,calls_for_tree[i,]))+"\n")
    else:
        if indels:
            subprocess.run(["rm -fr tree_counting/indels ; mkdir -p tree_counting/indels" ],shell=True)
        else: 
            subprocess.run(["rm -fr tree_counting/snps ; mkdir -p tree_counting/snps" ],shell=True)
    # build table    
    if indels:
        with open("tree_counting/indels/for_tree_labeling.csv",'w') as csv_file:
            csv_file.write(",".join(np.append(np.array(['chr','pos']),treeSampleNamesLong))+"\n") # write header
            for i in range(p_chr_table.shape[0]):
                csv_file.write(",".join(np.append( np.array([str(p_chr_table[i,0]),str(p_chr_table[i,1])]) ,calls_for_tree[i,]))+"\n")
    else:
        with open("tree_counting/snps/for_tree_labeling.csv",'w') as csv_file:
            csv_file.write(",".join(np.append(np.array(['chr','pos']),treeSampleNamesLong))+"\n") # write header
            for i in range(p_chr_table.shape[0]):
                csv_file.write(",".join(np.append( np.array([str(p_chr_table[i,0]),str(p_chr_table[i,1])]) ,calls_for_tree[i,]))+"\n")


# SNP counting
def count_number_mutational_events(ancestral_reconstruction_labelled_tree,ancestral_reconstruction_fasta, skip_preterminals=True, track_nodes=False):
    """
    Counts the number of mutational events at all positions.
    Uses a treetime output ancestral reconstruction tree (with labelled internal nodes) and calls for all the goodpos that went into building the tree and ancestral reconstruction. 
    These outputs are from treetime ancestral, see function create_ancestral_reconstruction.
    
    Homoplasic locations are all cells with value >1

    Options to count or skip terminal node transitions, (eg only count internal nodes).

    Default output names from treetime:
    annotated_tree.nexus == ancestral_reconstruction_labelled_tree
    ancestral_sequences.fasta == ancestral_reconstruction_fasta

    Example call:
    count_number_mutational_events(f'{analysis_params_output_name}_ancestral_reconstruction/annotated_tree.nexus',f'{analysis_params_output_name}_ancestral_reconstruction/ancestral_sequences.fasta'"""
    # parse fasta
    calls_pos={}
    for record in SeqIO.parse(ancestral_reconstruction_fasta, 'fasta'):
        calls_pos[record.id] = str(record.seq)
    # parse tree structure into adjacency graph to traverse
    parsed_tree = parse_tree(ancestral_reconstruction_labelled_tree) 
    root=parsed_tree.common_ancestor(parsed_tree.get_terminals())
    net = Phylo.to_networkx(parsed_tree)
    tree_as_dict_of_lists=networkx.to_dict_of_lists(net)
    # start tree traversal, checking if each internal node has same call as parent, if not iterate the val for that base +1
    transitions_goodpos=np.zeros(len(str(record.seq)))
    to_visit=[(x, root) for x in tree_as_dict_of_lists[root]]
    visited=set()
    visited.add(root)
    nodes_output={i:[] for i in range(len(transitions_goodpos))}
    while len(to_visit)>0:
        currently_processing=to_visit[0]
        parent=currently_processing[1]
        current_node=currently_processing[0]
        visited.add(current_node)
        is_preterminal=len(tree_as_dict_of_lists[current_node])==1
        if skip_preterminals and is_preterminal:
            pass
        else:
            for index in range(len(transitions_goodpos)):
                if calls_pos[parent.name][index] != calls_pos[current_node.name][index] and calls_pos[current_node.name][index]!= 'N' and calls_pos[parent.name][index]!= 'N':
                    transitions_goodpos[index]+=1
                    if track_nodes:
                        nodes_output[index].append((current_node.name,calls_pos[current_node.name][index]))
            for children in tree_as_dict_of_lists[current_node]:
                if children not in visited:
                    to_visit.append((children,current_node))
        to_visit=to_visit[1:]
    if track_nodes:
        return transitions_goodpos,nodes_output
    return transitions_goodpos

###################################################################################################
# # # # # # # # # # # # # # # # Geographic-related analyses # # # # # # # # # # # # # # # # # # # #
###################################################################################################
def measure_dist_from_coords(lat1, lon1, lat2, lon2):
    """Generally used geo measurement function
        Returns value in meters
    """
    radius = 6378.137 # Radius of earth in KM
    dLat = lat2 * math.pi / 180 - lat1 * math.pi / 180
    dLon = lon2 * math.pi / 180 - lon1 * math.pi / 180
    a = (math.sin(dLat/2) * math.sin(dLat/2) + math.cos(lat1 * math.pi / 180) * math.cos(lat2 * math.pi / 180) * math.sin(dLon/2) * math.sin(dLon/2))
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
    return d * 1000

