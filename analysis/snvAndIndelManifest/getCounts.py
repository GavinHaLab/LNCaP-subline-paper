#!/usr/bin/env python3

import sys, os, pandas as pd, numpy as np, re

## so, need to read in all samples vcf file. For each sample, obtain counts for total number of snvs, then counts per category: 3'UTR, 5'UTR, Intergenic_region, Intron, Missense, Start_lost, Stop_gained, Synonymous, Pathogenic

## read in vcf file
vcf = pd.read_csv(sys.argv[1], sep='\t', comment='#', header=0)

# need each of these counts per sample, so easiest is to loop through each sample, and then count the number of each type of mutation and fill in all categories in a final output df
outDF = pd.DataFrame(columns=['Sample', 'Total', '3UTR', '5UTR', 'Intergenic_region', 'Intron', 'Missense', 'Start_lost', 'Stop_gained', 'Synonymous', 'Pathogenic', 'snvCoding', 'indelCoding', 'Frameshift'])

# get all sample names
samples = vcf['Sample'].unique()
print(samples)

for sample in samples:
    # print(sample)
    sampleDF = vcf[vcf['Sample'] == sample]
    # print(sampleDF)
    total = sampleDF.shape[0]
    # print(total)

    # now get count for values in Func.refGene
    allFuncs = sampleDF['Func.refGene'].unique()    
    
    # now count each type of mutation
    muts = sampleDF['Mutation_Type'].unique()
    # print(muts)
    counts = dict()
    for mut in muts:
        counts[mut] = sampleDF[sampleDF['Mutation_Type'] == mut].shape[0]

    for func in allFuncs:
        counts[func] = sampleDF[sampleDF['Func.refGene'] == func].shape[0]

    ## to assess pathogenicity, we need to see if InterVar or ClinVar has a value of Pathogenic... do NOT double count
    ## need to make one column that has a combined pathogenicity value across InterVar, ClinVar, and Cosmic70. Make it T/F, if any of these columns has a value of Pathogenic or likely pathogenic, then this column should be True
    sampleDF['Pathogenic'] = (sampleDF['InterVar'] == 'Pathogenic') | (sampleDF['InterVar'] == 'Likely pathogenic') | (sampleDF['ClinVar'] == 'Pathogenic') | (sampleDF['ClinVar'] == 'Likely pathogenic') | (sampleDF['Cosmic70'] != '.')
    pathogenic = sampleDF[sampleDF['Pathogenic'] == True].shape[0]
    counts['Pathogenic'] = pathogenic
    counts['snvCoding'] = sampleDF[(sampleDF['Func.refGene']==('exonic')) | (sampleDF['Func.refGene']==('exonic;splicing')) & (sampleDF['Type'] == 'SNV')].shape[0] ## being a bit more specific here, not just looking for only exonic substring
    counts['indelCoding'] = sampleDF[(sampleDF['Func.refGene']==('exonic')) | (sampleDF['Func.refGene']==('exonic;splicing')) & (sampleDF['Type'] == 'INDEL')].shape[0]

    ## need to get frameshifts, insertions, deletions, etc. from AA_Change...
    counts['Frameshift'] = sampleDF[sampleDF['AA_Change'].str.contains('fs')].shape[0]
    counts['Insertion'] = sampleDF[sampleDF['AA_Change'].str.contains('ins')].shape[0]
    counts['Deletion'] = sampleDF[sampleDF['AA_Change'].str.contains('del')].shape[0]

    counts['nonSynonymous'] = sampleDF[sampleDF['Mutation_Type'].str.contains('nonsynonymous')].shape[0]
    counts['splicing'] = sampleDF[sampleDF['Func.refGene'].str.contains('splicing')].shape[0]
    counts['ncRNA'] = sampleDF[sampleDF['Func.refGene'].str.contains('ncRNA')].shape[0]
    counts['ncRNA_exonic'] = sampleDF[sampleDF['Func.refGene'].str.contains('ncRNA_exonic')].shape[0]

    # now fill in the final output df
    outDF = outDF.append({'Sample': sample, 'Total': total, '3UTR': counts.get('UTR3', 0), '5UTR': counts.get('UTR5', 0), 'Intergenic_region': counts.get('intergenic', 0), 'Intron': counts.get('intronic', 0), 'Missense': counts.get('nonsynonymous SNV', 0), 'Start_lost': counts.get('startloss', 0), 'Stop_gained': counts.get('stopgain', 0), 'Synonymous': counts.get('synonymous SNV', 0), 'Pathogenic': counts.get('Pathogenic', 0), 'snvCoding': counts.get('snvCoding', 0), 'indelCoding': counts.get('indelCoding', 0), 'Frameshift': counts.get("Frameshift", 0), 'nonSynonymous': counts.get('nonSynonymous', 0), 'splicing': counts.get('splicing', 0), 'ncRNA': counts.get('ncRNA', 0), 'ncRNA_exonic': counts.get('ncRNA_exonic', 0)}, ignore_index=True)

print(outDF)
outDF.to_csv(sys.argv[2], index=False, sep='\t')

## would like to transpose this df so that each row is a sample and each column is a category
outDF = outDF.set_index('Sample').T
## reorder columns to be the following: LNCaP_FGC	LNCaP_clone_FGC	LNCaP_C4	LNCaP_C42	LNCaP_C42B	LNCaP_95	LNCaP_Abl	LNCaP_16D	LNCaP_42D	LNCaP_42F	LNCaP_shAR-pATK	LNCaP_APIPC	LNCaP-AR-907	LNCaP-AR-909
outDF = outDF[['LNCaP_FGC', 'LNCaP_clone_FGC', 'LNCaP_C4', 'LNCaP_C42', 'LNCaP_C42B', 'LNCaP_95', 'LNCaP_Abl', 'LNCaP_16D', 'LNCaP_42D', 'LNCaP_42F', 'LNCaP_shAR-pATK', 'LNCaP_APIPC', 'LNCaP-AR-907', 'LNCaP-AR-909', 'LNCaP_FGC_SRR7943697']]
outDF.to_csv(sys.argv[3], sep='\t')
