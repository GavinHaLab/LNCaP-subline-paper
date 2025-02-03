#!/usr/bin/env python3

## read in cnv calls, reformat and sort by sample after fixing sample names

import pandas as pd, numpy as np, os, sys

# inMat = sys.argv[1]
outMat = sys.argv[1]

sample_list=['LNCaP_16D',
    'LNCaP_42D',
    'LNCaP_42F',
    'LNCaP_95',
    'LNCaP_Abl',
    'LNCaP_APIPC',
    'LNCaP_C4',
    'LNCaP_C42',
    'LNCaP_C42B',
    'LNCaP_FGC_PRJNA361316',
    'LNCaP_FGC',
    'LNCaP_shAR',
    'LNCaP-AR-907',
    'LNCaP-AR-909', 
    'LNCaP_FGC_SRR7943697'
    ]

comut_CN = pd.read_csv('path/to/comut_CN_cf0.8.csv', sep=',') ## output from LNCaP-subline-paper/analysis/makeMatrixFromTitan/generate_comut_CN.R

comut_CN = comut_CN.drop(['Unnamed: 0'], axis=1)
## Re-name the samples to their un-titan forms
comut_CN['sample'] = comut_CN['sample'].replace(['LNCaP_16D_recalibrated_cluster3',
    'LNCaP_42D_recalibrated_cluster1',
    'LNCaP_42F_recalibrated_cluster1',
    'LNCaP_95_recalibrated_cluster2',
    'LNCaP_Abl_recalibrated_cluster2',
    'LNCaP_APIPC_recalibrated_cluster3',
    'LNCaP_C4_recalibrated_cluster3',
    'LNCaP_C42_recalibrated_cluster2',
    'LNCaP_C42B_recalibrated_cluster3',
    'LNCaP_FGC_PRJNA361316_recalibrated_cluster1',
    'LNCaP_FGC__recalibrated_cluster3',
    'LNCaP_shAR-pATK_recalibrated_cluster1',
    'LNCaP-AR-907_cluster1',
    'LNCaP-AR-909_cluster2',
    'LNCaP_FGC_SRR7943697_cluster2'],
sample_list)

## sort by sample then by gene... sample, category
comut_CN = comut_CN.sort_values(by=['sample', 'category'])

comut_CN.to_csv(outMat, sep='\t', index=False)