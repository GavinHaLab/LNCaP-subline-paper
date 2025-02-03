#!/usr/bin/env python3

import os, sys, pandas as pd, numpy as np, argparse

def getArgs():
    parser = argparse.ArgumentParser(description='Get counts of unique and shared variants in a concatenated vcf file')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input concatenated manifest file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output counts matrix')
    return parser.parse_args()

args = getArgs()
infile = args.input
o = args.output

## rewrite in order to improve functionality when taking in a file with all samples and subsetting.

## read in input, then use count columns to assess the number of variants shared by all, shared between two, and unique to each
samples = ['LNCaP_FGC_PRJNA361316', 'LNCaP-FGC', 'LNCaP_FGC_SRR7943697']

## basic goal: subset columns to Gene	Start-Stop	Alt Base	Cytoband	Exon:AminoAcidChange:ProteinChange	Type	[samples]	Grand Total	Gene_Start_Stop_AltBase_Cytoband_Exon:AminoAcidChange:ProteinChange

cols = ['Gene', 'Start-Stop', 'Alt Base', 'Cytoband', 'Exon:AminoAcidChange:ProteinChange', 'Type'] + samples + ['Grand Total', 'Gene_Start_Stop_AltBase_Cytoband_Exon:AminoAcidChange:ProteinChange']

## want only unique comparisons... should be 3
comparisons = []
for sample1 in samples:
    for sample2 in samples:
        if sample1 == sample2:
            continue
        if [sample2, sample1] in comparisons:
            continue
        comparisons.append([sample1, sample2])

print(comparisons)
# exit()

## now need to map naming of comparisons (s1-VERSUS-s2) to the appropriate columns in the output matrix
colMap = dict()
for i in range(len(comparisons)):
    colMap[f'{comparisons[i][0]}-VERSUS-{comparisons[i][1]}'] = comparisons[i]

print(colMap)

inDf = pd.read_csv(infile, sep='\t', header=0)

## subset to only the columns we need
inDf = inDf[cols]

## recompute Grand Total by summing sample columns... fill NA with 0
inDf['Grand Total'] = inDf[samples].sum(axis=1)

## drop any rows where Grand Total is 0
inDf = inDf[inDf['Grand Total'] > 0]


inDf.to_csv('input.reformatted.txt', sep='\t', header=True, index=False)

if len(samples) > 2:
    ## print number of rows in input
    print(inDf.shape[0])
    unique = inDf[inDf['Grand Total'] == 1]
    shared = inDf[inDf['Grand Total'] == 2]
    allShared = inDf[inDf['Grand Total'] == 3].shape[0]

    ## output matrix will have rows: samples (one per sample for unique), allShared, colMap.keys()
    outDf = pd.DataFrame(index=samples + ['allShared'] + list(colMap.keys()), columns=['count'])

    for s in samples:
        outDf.loc[s, 'count'] = unique[unique[s] == 1].shape[0]

    outDf.loc['allShared', 'count'] = allShared

    for c in colMap.keys():
        print(c)
        print(colMap[c])
        s1 = colMap[c][0]
        s2 = colMap[c][1]
        outDf.loc[c, 'count'] = shared[(shared[s1] == 1) & (shared[s2] == 1)].shape[0]

    outDf.to_csv(o, sep='\t', header=True, index=True)

else:
    ## print number of rows in input
    ## replace NA with 0 in each sample column
    inDf = inDf.fillna(0)

    print(inDf.shape[0])
    unique = inDf[inDf['Grand Total'] == 1]
    shared = inDf[inDf['Grand Total'] == 2]

    ## output matrix will have rows: samples (one per sample for unique), allShared, colMap.keys()
    outDf = pd.DataFrame(index=samples + list(colMap.keys()), columns=['count'])
    print(outDf)
    print(colMap.keys())
    print(unique)
    for s in samples:
        # print(unique[unique[s ==1]])
        unique
        outDf.loc[s, 'count'] = unique[unique[s] == 1].shape[0]

    ## here we just need to get the shared counts since there is only one comparison
    outDf[colMap.keys()[0]] = shared.shape[0]

    outDf.to_csv(o, sep='\t', header=True, index=True)