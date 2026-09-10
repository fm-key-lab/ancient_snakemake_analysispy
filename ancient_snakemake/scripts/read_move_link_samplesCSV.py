# coding: utf-8

import os
from pathlib import Path
import sys
import glob
import subprocess
import pandas as pd
import numpy as np


def read_samplesCSV(spls):
    header_check = ['Path', 'Sample', 'Reference', 'Callindels', 'Outgroup', 'Ancient']
    parsed_samples=pd.read_csv(spls,sep=',',header=0,dtype=str)
    if list(parsed_samples.columns)!=header_check:
        raise TypeError(f'Header incorrect, should follow format : {",".join(header_check)}')
    if parsed_samples['Sample'].duplicated().any():
        raise ValueError('Sample names must be unique')
    for column in ['Callindels', 'Outgroup', 'Ancient']:
        if not parsed_samples[column].isin(['0', '1']).all():
            raise ValueError(f'{column} must contain only boolean values encoded as 0 or 1')
    numpy_parsed=parsed_samples.to_numpy()
    paths,samples,references,call_indels,outgroup,ancient=numpy_parsed.T
    # confirm path exists on all paths
    collector=[]
    for path in paths:
        if not os.path.isfile(path):
            collector.append(path)
    if len(collector)>0:
        raise ValueError(f'Paths not found for following paths {collector}')
    makelink_ancient(paths,samples,references)
    generate_freebayes_input(samples, references, call_indels, outgroup)
    return [paths,samples,references,call_indels,outgroup,ancient]

def parse_multi_genome_smpls(SAMPLE_ls,REF_Genome_ls):
    ## Expand lists if multiple genomes are used within a sample
    REF_Genome_ext_ls = []
    SAMPLE_ext_ls = []
    for sampleID, refgenomes in zip(SAMPLE_ls, REF_Genome_ls):
        for refgenome in refgenomes.split(" "):
            REF_Genome_ext_ls.append(refgenome)
            SAMPLE_ext_ls.append(sampleID)
    return [REF_Genome_ext_ls, SAMPLE_ext_ls]

def makelink_ancient(paths,samples,references):
    for bam,sample,reference in zip(paths,samples,references):
        os.makedirs(f'data/{reference}/{sample}', exist_ok=True)
        subprocess.run(f'ln -s -T {bam} data/{reference}/{sample}/{sample}.bam || echo data/{reference}/{sample}/{sample}.bam link path already exists, skipping', shell=True)
        subprocess.run(f'ln -s -T {bam}.bai data/{reference}/{sample}/{sample}.bam.bai || echo data/{reference}/{sample}/{sample}.bam.bai link path already exists, skipping', shell=True)

def get_bams(SAMPLE_ls, REF_Genome_ls):
    ## note: multiple ref genomes + varied ingroup/outgroup identity might be edgecase, if a sample is ingrp for one ref, outgrp for the other
    bams_ls = {}
    for sampleID,refgenomes in zip(SAMPLE_ls,REF_Genome_ls):
        for refgenome in refgenomes.split(" "):
            if refgenome not in bams_ls:
                bams_ls[refgenome] = [sampleID]
            else: 
                bams_ls[refgenome].append(sampleID)
    return bams_ls

def generate_freebayes_input(SAMPLE_ls, REF_Genome_ls, CALLINDELS_ls, OUTGROUP_ls):
    ## note: multiple ref genomes + varied ingroup/outgroup identity might be edgecase, if a sample is ingrp for one ref, outgrp for the other
    os.makedirs('0-freebayes_input/', exist_ok=True)
    non_outgroup_sample_ls = {x:[] for x in np.unique(REF_Genome_ls)}
    for sampleID,refgenomes,outgroup_bool,call_indels_bool in zip(SAMPLE_ls,REF_Genome_ls,OUTGROUP_ls,CALLINDELS_ls):
        for refgenome in refgenomes.split(" "):
            if int(outgroup_bool) == 0 and int(call_indels_bool) == 0:
                non_outgroup_sample_ls[refgenome].append(sampleID)
    # check if file exists, overwriting would restart the processing of indels
    for refgenome in non_outgroup_sample_ls:
        this_reference_output=f'0-freebayes_input/ref_{refgenome}_non_outgroup_bams.txt'
        if os.path.isfile(this_reference_output):
            with open(this_reference_output,'r') as f:
                contents=[l.strip() for l in f.readlines()]
            if non_outgroup_sample_ls[refgenome]==contents:
                break
        with open(f'0-freebayes_input/ref_{refgenome}_non_outgroup_bams.txt', "w") as f:
            for bam in non_outgroup_sample_ls[refgenome]:
                f.write(f"{bam}\n")
