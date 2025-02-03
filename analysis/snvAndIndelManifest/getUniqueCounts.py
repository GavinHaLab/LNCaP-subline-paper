#!/usr/bin/env python3

import pandas as pd, numpy as np, os, sys, re
import argparse

## need to read in concatenated vcf file, get count of variants unique to each sample and variants shared by all samples

def getArgs():
    parser = argparse.ArgumentParser(description='Get counts of unique and shared variants in a concatenated vcf file')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input concatenated vcf file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
    return parser.parse_args()

def main():
    ## get arguments
    args = getArgs()
    vcf = args.input
    out = args.output

    ## read in vcf file
    df = pd.read_csv(vcf, sep='\t', comment='#', header=0)

    ## how do I want to do this? first, split on caller -- mutect and svaba, which will separate snvs from indels... but how to get unique variants per sample?

    ## first, split on caller
    # mutect = df[df['Source'] == 'mutect']
    # svaba = df[df['Source'] == 'svaba']

    mutect = df[df['Type'] == 'SNV']
    svaba = df[df['Type'] == 'INDEL']


    ## changed... no longer using svaba. just using mutect, so split on snv/indel in Type column

    ## another question: how would the count of variant shared by all change per line?? different values in the table seem tos how that it could, but not sure how that's defined...
    ## confirmed -- shared by all is just the set of variants found in all samples, so it's the same for each line. one value per line, not per sample.

    ## get counts of unique variants per sample
    mutect['varID'] = mutect['CHR'].astype(str) + '_' + mutect['Start'].astype(str) + '_' + mutect['REF'] + '_' + mutect['ALT']
    svaba['varID'] = svaba['CHR'].astype(str) + '_' + svaba['Start'].astype(str) + '_' + svaba['REF'] + '_' + svaba['ALT']

    ## for each sample, get the count of total variants
    mutect_total = mutect.groupby('Sample')['varID'].count().reset_index()
    svaba_total = svaba.groupby('Sample')['varID'].count().reset_index()
    print(mutect_total)
    print(svaba_total)

    mutect_vars = mutect.groupby('Sample')['varID'].apply(set)
    svaba_vars = svaba.groupby('Sample')['varID'].apply(set)

    ## print out counts per sample from mutect_vars and svaba_vars... quick QC
    for sample in mutect_vars.index:
        print(f'{sample}: {len(mutect_vars[sample])}')
    for sample in svaba_vars.index:
        print(f'{sample}: {len(svaba_vars[sample])}')

    unique_counts_mutect = pd.Series({
        sample: len(variants - set.union(*mutect_vars.drop(sample, errors='ignore')))
        for sample, variants in mutect_vars.items()
    })

    unique_counts_svaba = pd.Series({
        sample: len(variants - set.union(*svaba_vars.drop(sample, errors='ignore')))
        for sample, variants in svaba_vars.items()
    })

    shared_mutect = set.intersection(*mutect_vars)
    shared_svaba = set.intersection(*svaba_vars)

    countSharedMutect = len(shared_mutect)
    countSharedSvaba = len(shared_svaba)

    print(unique_counts_mutect)
    print(unique_counts_svaba)
    print(countSharedMutect)
    print(countSharedSvaba)

    ## also get exonic-only counts, so subset to only those annotated as exonic
    # mutect_exonic = mutect[mutect['Func.refGene'] == 'exonic']
    # svaba_exonic = svaba[svaba['Func.refGene'] == 'exonic']
    mutect_exonic = mutect[mutect['Func.refGene'].isin(['exonic', 'exonic;splicing'])]
    svaba_exonic = svaba[svaba['Func.refGene'].isin(['exonic', 'exonic;splicing'])]


    mutect_exonic_vars = mutect_exonic.groupby('Sample')['varID'].apply(set)
    svaba_exonic_vars = svaba_exonic.groupby('Sample')['varID'].apply(set)

    unique_counts_mutect_exonic = pd.Series({
        sample: len(variants - set.union(*mutect_exonic_vars.drop(sample, errors='ignore')))
        for sample, variants in mutect_exonic_vars.items()
    })

    unique_counts_svaba_exonic = pd.Series({
        sample: len(variants - set.union(*svaba_exonic_vars.drop(sample, errors='ignore')))
        for sample, variants in svaba_exonic_vars.items()
    })

    shared_mutect_exonic = set.intersection(*mutect_exonic_vars)
    shared_svaba_exonic = set.intersection(*svaba_exonic_vars)

    countSharedMutectExonic = len(shared_mutect_exonic)
    countSharedSvabaExonic = len(shared_svaba_exonic)

    print(unique_counts_mutect_exonic)
    print(unique_counts_svaba_exonic)
    print(countSharedMutectExonic)
    print(countSharedSvabaExonic)

    ## write to output file
    with open(out, 'w') as f:
        f.write('Sample\tMutect_Total\tSvaba_Total\tMutect_Unique\tSvaba_Unique\tShared_Mutect\tShared_Svaba\tMutect_Exonic\tSvaba_Exonic\tShared_Mutect_Exonic\tShared_Svaba_Exonic\n')
        for sample in unique_counts_mutect.index:
            f.write(f'{sample}\t{mutect_total[mutect_total["Sample"] == sample]["varID"].values[0]}\t{svaba_total[svaba_total["Sample"] == sample]["varID"].values[0]}\t{unique_counts_mutect[sample]}\t{unique_counts_svaba[sample]}\t{countSharedMutect}\t{countSharedSvaba}\t{unique_counts_mutect_exonic[sample]}\t{unique_counts_svaba_exonic[sample]}\t{countSharedMutectExonic}\t{countSharedSvabaExonic}\n')
    ## read in the output as matrix, transpose so samples are columns, and then write to a new file
    outDF = pd.read_csv(out, sep='\t', header=0)
    print(out)
    out = out.replace('.txt', '_transposed.txt')
    print(out)
    outDF = outDF.set_index('Sample').T
    outDF = outDF[['LNCaP_FGC', 'LNCaP_clone_FGC', 'LNCaP_C4', 'LNCaP_C42', 'LNCaP_C42B', 'LNCaP_95', 'LNCaP_Abl', 'LNCaP_16D', 'LNCaP_42D', 'LNCaP_42F', 'LNCaP_shAR-pATK', 'LNCaP_APIPC', 'LNCaP-AR-907', 'LNCaP-AR-909', 'LNCaP_FGC_SRR7943697']]

    outDF.to_csv(out, sep='\t', index=True)

if __name__ == '__main__':
    main()