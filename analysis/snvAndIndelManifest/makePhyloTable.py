#!/usr/bin/env python3

import pandas as pd, numpy as np, os, sys, re
import argparse
## new approach: read in manifest file, get unique values for var_id (CHR_Start_REF_ALT), then determine which samples have each of these variants...
## end result will be a table with the following columns: Gene,Start-Stop,Alt Base,Cytoband,Exon:AminoAcidChange:ProteinChange,Type,LNCaP_16D,LNCaP_42D,LNCaP_42F,LNCaP_95,LNCaP_Abl,LNCaP_APIPC,LNCaP_C4,LNCaP_C42,LNCaP_C42B,LNCaP_dbGaP,LNCaP_FGC,LNCaP_shAR-pATK,Grand Total,Gene_Start_Stop_AltBase_Cytoband_Exon:AminoAcidChange:ProteinChange[Cosmic70][InterVar][ClinVar]

## columns in input are:
# Sample CHR Start End REF ALT Func.refGene Gene Mutation_Type Cosmic70 InterVar ClinVar cytoBand AD_REF AD_ALT DP VAF Source Type AA_Change
## so, read in the df containing all of these values, get unique varID values, make output df

def getArgs():
    parser = argparse.ArgumentParser(description='Get counts of unique and shared variants in a concatenated vcf file')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input concatenated manifest file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
    return parser.parse_args()

args = getArgs()
i = args.input
o = args.output
inDf = pd.read_csv(i, sep='\t', header=0)

## get unique varID values
inDf['varID'] = inDf['CHR'].astype(str) + '_' + inDf['Start'].astype(str) + '_' + inDf['REF'] + '_' + inDf['ALT']
# samplesAllowed = ['LNCaP_FGC', 'LNCaP_clone_FGC', 'LNCaP_FGC_SRR7943697']
# inDf = inDf[inDf['Sample'].isin(samplesAllowed)]
## print sample names
print(inDf['Sample'].unique())
inDf['varID'] = inDf['varID'].astype(str)


## print unique counts for Mutation_Type
inDf = inDf[inDf['Func.refGene'].isin(['exonic', 'exonic;splicing'])]
inDf.loc[inDf['Mutation_Type'] == '.', 'Mutation_Type'] = inDf.loc[inDf['Mutation_Type'] == '.', 'Func.refGene']
# print(inDf['varID'].dtype)
uniqueVars = inDf['varID'].unique()

varFuncs = dict()
inDf['Start-Stop'] = inDf['Start'].astype(str) + '-' + inDf['End'].astype(str)
inDf['Exon:AminoAcidChange:ProteinChange'] = inDf['AA_Change']

# Select relevant columns
relevant_columns = {
    'Gene': 'Gene',
    'Start-Stop': 'Start-Stop',
    'ALT': 'Alt Base',
    'cytoBand': 'Cytoband',
    'Exon:AminoAcidChange:ProteinChange': 'Exon:AminoAcidChange:ProteinChange',
    'Mutation_Type': 'Type'
}

# Aggregate data by 'varID' and define an aggregation function
aggregated_df = (
    inDf.groupby('varID')
    .agg({col: 'first' for col in relevant_columns.keys()})
    .rename(columns=relevant_columns)
)

# Convert the aggregated DataFrame to dictionary
varFuncs = aggregated_df.to_dict(orient='index')

print('done making dict')

cols = ['varID'] + list(inDf['Sample'].unique()) + ['Grand Total']
print(cols)
outDf = pd.DataFrame(columns=cols, index=uniqueVars)
## set varID values to be the index
outDf['varID'] = uniqueVars
i = 0
for var in uniqueVars:
    i += 1
    if i % 10000 == 0:
        print(i, ' variants processed so far')
    samples = inDf[inDf['varID'] == var]['Sample']
    for sample in samples:
        outDf.loc[outDf['varID'] == var, sample] = 1
        ### if other values not present (aa_change, etc.), add them here, but not grand total. Will have to do this later.
    outDf.loc[outDf['varID'] == var, 'Grand Total'] = samples.shape[0]
    outDf.loc[outDf['varID'] == var, 'Gene'] = varFuncs[var]['Gene']
    outDf.loc[outDf['varID'] == var, 'Start-Stop'] = varFuncs[var]['Start-Stop']
    outDf.loc[outDf['varID'] == var, 'Alt Base'] = varFuncs[var]['Alt Base']
    outDf.loc[outDf['varID'] == var, 'Cytoband'] = varFuncs[var]['Cytoband']
    outDf.loc[outDf['varID'] == var, 'Exon:AminoAcidChange:ProteinChange'] = varFuncs[var]['Exon:AminoAcidChange:ProteinChange']
    outDf.loc[outDf['varID'] == var, 'Type'] = varFuncs[var]['Type']
    outDf.loc[outDf['varID'] == var, 'Gene_Start_Stop_AltBase_Cytoband_Exon:AminoAcidChange:ProteinChange'] = varFuncs[var]['Gene'] + '_' + varFuncs[var]['Start-Stop'] + '_' + varFuncs[var]['Alt Base'] + '_' + varFuncs[var]['Cytoband'] + '_' + varFuncs[var]['Exon:AminoAcidChange:ProteinChange']

outDf = outDf[['Gene', 'Start-Stop', 'Alt Base', 'Cytoband', 'Exon:AminoAcidChange:ProteinChange', 'Type'] + list(inDf['Sample'].unique()) + ['Grand Total', 'Gene_Start_Stop_AltBase_Cytoband_Exon:AminoAcidChange:ProteinChange']]

# nameDict = {'LNCaP_FGC': 'LNCaP-FGC', 'LNCaP_clone_FGC': 'LLNCaP_FGC_PRJNA361316', 'LNCaP_FGC_SRR7943697': 'LNCaP_FGC_SRR7943697'}

outDf.columns = [nameDict.get(x, x) for x in outDf.columns]

## sort by gene, A->Z
outDf = outDf.sort_values(by='Gene')


outDf.to_csv(o, sep='\t', index=False)

