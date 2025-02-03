#!/usr/bin/env python3

import pandas as pd, numpy as np, os, sys, argparse

def getArgs():
    parser = argparse.ArgumentParser(description = 'Get mutation frequency from a mutation matrix')
    parser.add_argument('-im', '--input_matrix', type = str, required = True, help = 'Input mutation matrix')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output mutation frequency file')

    args = parser.parse_args()
    return args

args = getArgs()

# Read in mutation matrix
mut_matrix = pd.read_csv(args.input_matrix, sep = '\t', header=0)

## get unique genes, genes in category column
genes = mut_matrix['category'].unique()

## get count of mutations per gene... subset first to only rows without 'No Mutation' and 'No CNA'... then keep only frameshift, nonsynonymous, and stopgain