#!/usr/bin/env python3

import pandas as pd, numpy as np, os, sys, argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Create input for comut from snv and indel manifest')
    parser.add_argument('-m', '--manifest', help='Manifest file with snv and indel calls', required=True)
    parser.add_argument('-g', '--genes', help='List of genes to include', required=True)
    parser.add_argument('-o', '--output', help='Output file name', required=True)
    return parser.parse_args()

args = parse_args()
m = args.manifest
o = args.output
g = args.genes

genesToKeep = [
    'AR', 'TP53', 'FOXA1', 'ZBTB16', 'NCOR1', 'NCOR2', 'PIK3CA', 'PIK3CB', 'PIK3R1', 'PTEN', 'AKT1', 'BRAF', 'APC', 'CTNNB1', 'RAAF1', 'RNF43', 'RSPO2', 'ZNRF3', 
    'BRCA2', 'ATM', 'BRCA1', 'CDK12', 'MLH1', 'MSH2', 'MSH6', 'RB1', 'CDKN1B', 'CDKN2A', 'CCND1', 'KMT2C', 'KMT2D', 'KDM6A', 'CHD1', 'SPOP', 'MED12', 'ZFHX3', 'ERF',
    'GNAS'
]

with open(g) as f:
    genes = []
    for line in f:
        gene = line.strip().split('\t')[3]
        if gene in genesToKeep:
            # genes.append(gene)
            print('Gene in list of genes to keep: {}'.format(gene))
        else:
            print('Gene not in list of genes to keep: {}'.format(gene))
            genes.append(gene)

print(genes)

df = pd.read_csv(m, sep='\t')
## subset to exonic
df = df[df['Func.refGene'].str.contains('exonic')]
# cosmic70 != ., clinvar != . or likely benign or benign or uncertain significance, intervar != . or benign or likely benign or uncertain significance
df_patho = df[(df['Cosmic70'] != '.') | (df['ClinVar'].isin(['Pathogenic', 'Likely pathogenic', 'Pathogenic/', 'Likely_pathogenic'])) | (df['InterVar'].isin(['Pathogenic', 'Likely pathogenic']))]
outIntermediate = o.split('.')[0] + '_intermediate.txt'
# df_path = df[df['Gene'].isin(genes)]
df_patho = df_patho[df_patho['Gene'].isin(genes)]
df_patho.to_csv(outIntermediate, sep='\t', index=False)
# df = df[df['Cosmic70'] != '.' | df['ClinVar'].isin(['Pathogenic', 'Likely pathogenic', 'Pathogenic/', 'Likely_pathogenic']) | df['InterVar'].isin(['Pathogenic', 'Likely pathogenic'])] ## this isn't working...
df = df_patho[['Sample', 'Gene', 'Mutation_Type']]
# df = df[df['Gene'].isin(genes)]



## rename columns for comut: sample,category,value
df.columns = ['sample', 'category', 'value']

df.to_csv(o, sep='\t', index=False)

## for each gene, get the number of samples with a mutation in that gene
freqOut = o.split('.')[0] + '_freq.txt'

with open(freqOut, 'w') as f:
    for gene in genes:
        # numSamples = df[df['category'] == gene].shape[0]
        ## need to figure out the number of samples with a mutation in the gene
        numSamples = df[df['category'] == gene]['sample'].nunique()
        f.write('{}\t{}\n'.format(gene, numSamples))
# for gene in genes:
    


'''
python3 makeComutInput.py -m /path/to/snvAndIndel/manifest.txt -g /path/to/GeneList_for_comut.txt -o /path/to/output.txt 
'''