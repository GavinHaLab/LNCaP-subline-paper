#!/usr/bin/env python3

import os, sys, pandas as pd, numpy as np, argparse
import multiprocessing as mp
def getArgs():
    parser = argparse.ArgumentParser(description='Get counts of unique and shared variants in a concatenated vcf file')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input mutect2 run dir')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output vcf dir to be used in SigProfiler')
    parser.add_argument('-n', '--numCpu', type=int, required=False, default=mp.cpu_count(), help='Number of CPUs to use')
    return parser.parse_args()

args = getArgs()
i = args.input
o = args.output

## read in input per sample, filter to exclude germline/SNP events, output to properly formatted vcf files (one per sample)

## read in input
samples = os.listdir(i)

# for s in samples:
def getSampleInfo(i, s):
    vcf = os.path.join(i, s, 'pass_snvs.annovar.hg38_multianno.txt')
    if not os.path.exists(vcf):
        print(f'No file found for {s}')
        # continue
    inDf = pd.read_csv(vcf, sep='\t', header=0)
    ## filter to exclude germline/SNP events
    threshold = 0.1
    samp_df.loc[samp_df['gnomAD_genome_ALL'] == '.', 'gnomAD_genome_ALL'] = '0'
    out =  samp_df[pd.to_numeric(samp_df['gnomAD_genome_ALL']) > threshold].index
    filtered_df_init =  samp_df.drop(out)
    #### Remove Common Variants - EXAC ####
    filtered_df_init.loc[filtered_df_init['ExAC_ALL'] == '.', 'ExAC_ALL'] = '0'
    out = filtered_df_init[pd.to_numeric(filtered_df_init['ExAC_ALL']) > threshold].index
    filtered_df = filtered_df_init.drop(out)
    ## output to properly formatted vcf files

    filtered_df = filtered_df[['Chr', 'Start', 'Ref', 'Alt']]
    ## add sample column, fill with sample name. should go between start and ref
    filtered_df.insert(2, 'Sample', s)

    ## write to output, without header
    outVcf = os.path.join(o, f'{s}.vcf')
    filtered_df.to_csv(outVcf, sep='\t', header=False, index=False)

    # print(f'Wrote {filtered_df.shape[0]} variants to {outVcf}')

## parallelize

pool = mp.Pool(args.numCpu)
pool.starmap(getSampleInfo, [(i, s) for s in samples])
pool.close()
pool.join()


print('Done')
