#!/usr/bin/env python3
"""
Created on Thu Jan 13 15:04:55 2022

@author: ad_lemur
"""

#  Version History
# Tami Lieberman 2012, modified 2016
# Given a vcf file with one file per line, grabs the FQ score for each
# position. Ignores lines corresponding to indels
# Arolyn, 2018.12.18: Modified to work with snakemake!!
# Arolyn, 2018.12.18: Modified to work with gzipped files

import os
import gzip
import math 
import numpy as np
from analysispy_module import genomestats


def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5) / multiplier


def vcf_to_quals_snakemake(sample_path_to_vcf, sample_path_to_quals, REF_GENOME_DIRECTORY):

    ## Input file
    fname_in_gz = [os.getcwd() + '/' + sample_path_to_vcf] ## store zipped file names

    [ChrStarts,GenomeLength,ScafNames] = genomestats(REF_GENOME_DIRECTORY);

    print(f'\n {GenomeLength} \n') ## print to stdout or file??

    quals = np.zeros(GenomeLength, dtype=int)

    i = 1
    position = 0

    for file in fname_in_gz:
        with gzip.open(file, 'rt') as fid:

            for line in fid:
                if line[0] != "#":

                    ## print evry 50000 lines a dot (check for running status)
                    if not i%50000:
                        print('.', end = '')

                    ## parsing line
                    l = line.split('\t')

                    ## get single column values
                    chrom = np.where(ScafNames == l[0])
                    position_on_chr = int(l[1])
                    position = ChrStarts[chrom] + position_on_chr

                    ## store nucleotides
                    alt = l[4]
                    ref = l[3]

                    # only store quality for simple calls (not indel, not ambigious)
                    if (alt != '') & (',' not in alt) & (len(alt) == len(ref)) & (len(ref) == 1):
                        ## find and parse quality score in INFO column
                        info_col = l[7].split(';')

                        ## find FQ entry in info column and store value after '=' sign
                        entrywithFQ = [item for item in info_col if item.startswith('FQ')]

                        ## check if FQ was found
                        if not entrywithFQ:
                            print('No FQ entry was found in' + file + ' in row ' + str(i) + '!')
                        else:
                            ## Extract FQ value from string
                            entrywithFQ_str = entrywithFQ[0]
                            fq = float(entrywithFQ_str[entrywithFQ_str.index('=')+1:])

                            ## if already a position with a stronger FQ here, don't include this. (More negative is stronger)
                            if fq < quals[position-1]: ##Note python is 0 based!
                                quals[position-1] = int(round_half_up(fq))

                    ## Needed for progess bar on stdout
                    i += 1

    ## save file
    outpath = os.getcwd() + '/' + sample_path_to_quals
    np.savez_compressed(outpath, quals=quals)

    ## print to stdout
    print(f'\n Saved: {outpath} \n')
