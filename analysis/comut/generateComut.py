#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator # this function sets the location of the minor tick mark
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import palettable
from pathlib import Path
import os
import io

from comut import comut

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
    # 'LNCaP_FGC_SRR7943697'
    ]

# comut_CN = pd.read_csv('/path/to/comutFormattedCN/comut_CN_cf0.8.csv', sep=',')

# comut_CN = comut_CN.drop(['Unnamed: 0'], axis=1)
# ## Re-name the samples to their un-titan forms
# comut_CN['sample'] = comut_CN['sample'].replace(['LNCaP_16D_recalibrated_cluster3',
#     'LNCaP_42D_recalibrated_cluster1',
#     'LNCaP_42F_recalibrated_cluster1',
#     'LNCaP_95_recalibrated_cluster2',
#     'LNCaP_Abl_recalibrated_cluster2',
#     'LNCaP_APIPC_recalibrated_cluster3',
#     'LNCaP_C4_recalibrated_cluster3',
#     'LNCaP_C42_recalibrated_cluster2',
#     'LNCaP_C42B_recalibrated_cluster3',
#     'LNCaP_FGC_PRJNA361316_recalibrated_cluster1',
#     'LNCaP_FGC__recalibrated_cluster3',
#     'LNCaP_shAR-pATK_recalibrated_cluster1',
#     'LNCaP-AR-907_cluster1',
#     'LNCaP-AR-909_cluster2',
#     'LNCaP_FGC_SRR7943697_cluster2'],
# sample_list)

# comut_CN_filt = comut_CN.copy()
# comut_CN_freq = comut_CN_filt.groupby(['category']).agg('count').reset_index()
# comut_freq = comut_CN_freq.drop(['value'], axis=1)
# fin_comut_freq = comut_freq.rename(columns={'sample': 'Mutated samples'})
# fin_comut_freq

# ## Load in results from R pipeline of processing mutations
# var_genes_path = '/path/to/SNVandIndelCalls/comutFormatted/LNCaP_comutSNVinputs.txt'
# var_freq_path = '/path/to/SNVandIndelCalls/comutFormatted/LNCaP_comutSNVinputs_freq.txt'

# mutation = pd.read_csv(var_genes_path, sep='\t')
# # events = comut_CN
# mut_freq_df = pd.read_csv(var_freq_path, sep='\t')

# # mutation['value'].unique()

# # ## Drop non-frameshift or synonymous values
# functionals = ['nonsynonymous SNV', 'frameshift insertion', 'frameshift deletion', 'stopgain', 'startloss', 'splicing']

# mutations = (mutation.loc[mutation['value'].isin(functionals)]).reset_index(drop=True) # this will subset the data to only include the functional mutations... already mostly done in preprocessing script

# # for replacing values
# mutations['value'] = mutations['value'].replace(['startloss', 'splicing', 'stopgain', 'frameshift insertion', 'frameshift deletion', '.'], ['nonsynonymous SNV', 'nonsynonymous SNV', 'nonsense', 'frameshift', 'frameshift', 'nonsynonymous SNV'])

# ## MERGE comut_CN with mutation events
# new_muts = pd.merge(mutations, comut_CN, how='outer')
# ## sort new_muts by sample column, then category
# new_muts = new_muts.sort_values(by=['sample', 'category']).reset_index(drop=True)


# ## Check unique mutations and SV calls to organize for filling in nonexistent values
# CN_values = comut_CN['value'].unique().tolist()
# mutation_values = mutations['value'].unique().tolist()
# print("CN Unique values: ", comut_CN['value'].unique())
# print()
# print("Mutation Unique values: ", mutations['value'].unique())
## append ZZZZ to the beginning of 8p, 8q, 13q, 16p, 16q for sorting purposes
# value_dict = {'8q': 'ZZZZ8q', '8p': 'ZZZZ8p', '13q': 'ZZZZ13q', '16q': 'ZZZZ16q', '16p': 'ZZZZ16p'}

# genesToKeep = [
#     'AR', 'TP53', 'FOXA1', 'ZBTB16', 'NCOR1', 'NCOR2', 'PIK3CA', 'PIK3CB', 'PIK3R1', 'PTEN', 'AKT1', 'BRAF', 'APC', 'CTNNB1', 'RAAF1', 'RNF43', 'RSPO2', 'ZNRF3', 
#     'BRCA2', 'ATM', 'BRCA1', 'CDK12', 'MLH1', 'MSH2', 'MSH6', 'RB1', 'CDKN1B', 'CDKN2A', 'CCND1', 'KMT2C', 'KMT2D', 'KDM6A', 'CHD1', 'SPOP', 'MED12', 'ZFHX3', 'ERF',
#     'GNAS'
# ]

# ## need to add in rows for sample + gene combos that have no mutations
# inDF = new_muts
# cnv_values = ['Deletion (2 copies)', 'Gain', 'Amplification', 'Shallow Deletion (1 copy)', 'Deep Deletion (3 copies)', 'Copy_Neutral_LOH', 'Deletion LOH', 'Shallow Deletion LOH', 'Deep Deletion_LOH', 'CN_Homozygous_Del']
# snv_values = ['nonsynonymous SNV', 'frameshift', 'nonsense']

# for sample in sample_list:
#     for gene in genesToKeep:
#         # if values are present, check if there is one for cnv and one for snv using cnv_values and snv_values. if not, add in a row with the sample, gene, and value as 'No Mutation' or 'No CNA'
#         if inDF[(inDF['sample'] == sample) & (inDF['category'] == gene)].empty: ## if there is no value for the sample and gene
#             inDF = inDF.append({'sample': sample, 'category': gene, 'value': 'No Mutation'}, ignore_index=True)
#             inDF = inDF.append({'sample': sample, 'category': gene, 'value': 'No CNA'}, ignore_index=True)
#         else:
#             snv = False
#             cnv = False
#             if len(inDF[(inDF['sample'] == sample) & (inDF['category'] == gene)]) == 1:
#                 if inDF[(inDF['sample'] == sample) & (inDF['category'] == gene)]['value'].values[0] in snv_values:
#                     inDF = inDF.append({'sample': sample, 'category': gene, 'value': 'No CNA'}, ignore_index=True)
#                 elif inDF[(inDF['sample'] == sample) & (inDF['category'] == gene)]['value'].values[0] in cnv_values:
#                     inDF = inDF.append({'sample': sample, 'category': gene, 'value': 'No Mutation'}, ignore_index=True)
#             else:
#                 for value in inDF[(inDF['sample'] == sample) & (inDF['category'] == gene)]['value'].values:
#                     if value in snv_values:
#                         snv = True
#                     elif value in cnv_values:
#                         cnv = True
#                 if snv == False:
#                     inDF = inDF.append({'sample': sample, 'category': gene, 'value': 'No Mutation'}, ignore_index=True)
#                 if cnv == False:
#                     inDF = inDF.append({'sample': sample, 'category': gene, 'value': 'No CNA'}, ignore_index=True)
#         # snv = False
#         # cnv = False

# inDF['category'] = inDF['category'].replace(value_dict)

# inDF = inDF.sort_values(by=['sample', 'category']).reset_index(drop=True)
## now remove the ZZZZ from the beginning of the values
# inDF['category'] = inDF['category'].replace(['ZZZZ8q', 'ZZZZ8p', 'ZZZZ13q', 'ZZZZ16q', 'ZZZZ16p'], ['8q', '8p', '13q', '16q', '16p'])

# new_muts = inDF

new_muts = pd.read_csv('comutInputFilePath.txt', sep='\t', header=0)
freqs = pd.read_csv('inputMutationFrequencies.txt', sep='\t', header=0)
tmb = pd.read_csv('inputTmbValues.txt', sep='\t', header=0)
## read in group file
group_df2 = pd.read_csv('sampleGroupings.txt', sep='\t', header=0)
print(group_df2)




order = [
    'AR', 'TP53', 'PTEN', 'ETS Fusions', 'FOXA1', 'ZBTB16', 'NCOR1', 'NCOR2', 'PIK3CA', 'PIK3CB', 'PIK3R1', 'AKT1', 'BRAF', 'APC', 'CTNNB1', 'RAAF1', 'RNF43', 'RSPO2', 'ZNRF3', 
    'BRCA2', 'ATM', 'BRCA1', 'CDK12', 'MLH1', 'MSH2', 'MSH6', 'RB1', 'CDKN1B', 'CDKN2A', 'CCND1', 'KMT2C', 'KMT2D', 'KDM6A', 'CHD1', 'SPOP', 'MED12', 'ZFHX3', 'ERF',
    'GNAS', '8p', '8q',  '13q', '16p', '16q'
]


## new addition for SRR removal
order = order[::-1]
mutation = new_muts

func = ['nonsynonymous SNV', 'nonsense', 'frameshift']
freq_df = freqs
# sample_list = sample
## drop SRR from mutation

mutationNoSRR = mutation[~mutation['sample'].str.contains('SRR')]
# mutation3_tmp = mutation3[mutation3['value'].isin(func)]
muts_tmp = mutationNoSRR[mutationNoSRR['value'].isin(func)]

for gene in freq_df['category']:
  ## don't need shape, need count of unique sample ids
  print(gene)
  freq_df.loc[freq_df['category'] == gene, 'Mutation_frequency'] = muts_tmp[muts_tmp['category'] == gene]['sample'].nunique()
  print(muts_tmp[muts_tmp['category'] == gene])



## will get back to this later
mutation_burden3 = tmb
mutation_burdenNoSRR = mutation_burden3[~mutation_burden3['sample'].str.contains('SRR')]
percentages = ((freq_df['Mutation_frequency']/len(sample_list)*100).round(1).astype(str) + '%')
mutFreq_mapping = {'Mutation_frequency': 'darkgrey'}
# CNFreq_mapping = {'CN_frequency': 'darkgrey'}
group_df2 = group_df2[~group_df2['sample'].str.contains('SRR')]
## ALTERATION Mapping
alt_mapping = {'nonsynonymous SNV': '#FF796C', 'frameshift': 'plum', 'nonsense': 'slategrey', 
'Gain': {'facecolor': 'darkred'}, 
'Amplification': {'facecolor': 'red'}, 
'Shallow Deletion (1 copy)': {'facecolor': 'lightblue'}, 
## Added for last few rows
'>5% Arm Deleted': {'facecolor': 'lightblue'},
'>5% Arm Gained': {'facecolor': 'lightsalmon'},
'Fusion': {'facecolor': 'limegreen'},
'Unknown': {'facecolor': 'lightgrey'},
## Added for last few rows
'Deletion (2 copies)': {'facecolor': 'blue'},
'Deep Deletion (3 copies)': {'facecolor': 'purple'},
'Copy_Neutral_LOH': {'facecolor': 'none', 'hatch': '///', 'edgecolor': 'black', 'linewidth': 0.3},
'Deep Deletion_LOH': {'facecolor': 'purple', 'hatch': '///', 'edgecolor': 'yellow', 'linewidth': 0.3},
'Deletion LOH': {'facecolor': 'blue', 'hatch': '///', 'edgecolor': 'yellow', 'linewidth': 0.3},
'Shallow Deletion LOH': {'facecolor': 'lightblue', 'hatch': '///'},
'CN_Homozygous_Del': {'facecolor': 'darkgrey', 'hatch': '///', 'edgecolor': 'yellow', 'linewidth': 0.3},
'Absent': {'facecolor': 'lightgrey', 'alpha': 0.2},
'No CNA': {'facecolor': 'lightgrey', 'alpha': 0.2, 'edgecolor': 'lightgrey', 'linewidth': 0.3},
'No Mutation': {'facecolor': 'lightgrey', 'alpha': 0.2, 'edgecolor': 'lightgrey', 'linewidth': 0.3}, 
}
# 'Absent': {'facecolor': 'grey', 'alpha': 0.2}}
## SV-CN Mapping colors
# event_mapping = {'DEL_LOH': 'green', 'NEUTRAL': 'lightgrey', 'Gain': 'red', 'Amplification': 'cyan'}
event_mapping = {'Amplification': {'facecoor': 'green', 'edgecolor': 'black', 'linewidth': 1}, 'Shallow_Deletion': {'facecolor': 'orange', 'edgecolor': 'black', 'linewidth': 1}, 'Deep_Deletion': 'red'}
indicator_kwargs = {'color': 'black', 'marker': 'o', 'linewidth': 1, 'markersize': 5}
## ADD BAR PLOTS
# bar_mapping = {'Nonsynonymous': 'purple', 'Synonymous': 'pink'}
bar_mapping = {'Both': 'darkgrey'}
bar_kwargs = {'width': 0.3, 'edgecolor': 'black'}
# Add percentages to bar plot
## will get back to this later
# percentages = ((freq_df['Mutation_frequency']/len(sample_list)*100).round(1).astype(str) + '%')
mut_order = ['nonsynonymous SNV', 'frameshift', 'nonsense', 'Multiple', 'Amplification', 'Gain', 'Shallow Deletion (1 copy)', 'Shallow Deletion LOH', 'Deletion (2 copies)', 'Deletion LOH', 'Deep Deletion (3 copies)', 'Deep Deletion_LOH', 'CN_Homozygous_Del', 'Copy_Neutral_LOH', 
             '>5% Arm Gained', '>5% Arm Deleted', 'Fusion', 'Unknown']

# cn_order = ['Amplification', 'Gain', 'Shallow Deletion (1 copy)', 'Deletion (2 copies)', 'Deep Deletion (3 copies)']

value_order = ['Shallow Deletion (1 copy)', 'Deletion (2 copies)', 'Deletion LOH', 'Shallow Deletion LOH', 'CN_Homozygous_Del', 'Copy_Neutral_LOH', 'Deep Deletion_LOH', 'No CNA']

priority = ['Amplification', 'Deletion LOH', 'Shallow Deletion LOH', 'Deep Deletion LOH', 'Shallow Deletion (1 copy)', 'Deletion (2 copies)', 'Copy Neutral LOH', 'CN_Homozygous_Del', 'No CNA'']

################ BEGIN COMUT ##################

basal_comut = comut.CoMut()

# Add Grouping
basal_comut.add_sample_indicators(group_df2, name= 'Line Relation', plot_kwargs = indicator_kwargs)

#add mutation
basal_comut.add_categorical_data(mutationNoSRR, name = 'Event Type', mapping = alt_mapping, category_order = freq_df['category'], value_order = value_order, priority = priority, borders=borders)
# basal_comut.add_categorical_data(mutation, name = 'Event Type', mapping = alt_mapping, category_order = order, value_order = value_order, priority = priority)

#add mut frequency bar
## will get back to this later
basal_comut.add_side_bar_data(freq_df, paired_name = 'Event Type', name = 'Mutation_frequency', position = 'right', mapping = mutFreq_mapping, xlabel = 'Mutation Prevalence', bar_kwargs = {'alpha': 0.5, 'height': 0.8})

# basal_comut.add_side_bar_data(freq_df, paired_name = 'Event Type', name = 'CN_frequency', position = 'right', mapping = CNFreq_mapping, xlabel = 'Alteration Prevalence', bar_kwargs = {'alpha': 0.5, 'height': 0.8})


#Add BAR PLOTS
## will get back to this later
basal_comut.add_bar_data(mutation_burdenNoSRR, name='Mutation Burden', mapping=bar_mapping, stacked=False, bar_kwargs=bar_kwargs, ylabel='Mutations/Mb')


# global matplotlib params
custom_rcParams = {
  'font.family': 'Arial',
  'font.size': 10, #?
  'axes.titlesize' : 18, #?
  'axes.labelsize': 8, #plot label size
  'legend.fontsize': 10, #legend size 
  'xtick.labelsize': 10, #side bar plot x-axis ; samples
  'ytick.labelsize': 10 #gene names and y-axis bar plots 
}

rcParams.update(custom_rcParams)
border_white = ['Absent']

## ORIGINAL SIZE
# fig = plt.figure(figsize = (9, 14))
# widths = [0.3, 4] #le  ft/right side bar

## NEW SIZE
fig = plt.figure(figsize = (8,15))
widths = [2, 0.75] #le  ft/right side bar


#basal_comut.plot_comut(figsize = (25, 15), heights = heights, widths = widths, wspace = wspace, x_padding = 0.04, y_padding = 0.04, tri_padding = 0.03, hspace = 0.1)
basal_comut.plot_comut(x_padding = 0.03, y_padding = 0.03, widths = widths, tri_padding = 0.02,  subplot_hspace = 0.04, wspace = 0.04, hspace = 0.005, fig = fig)

# add a vertical line at -log(Q) = 1 and a text box to the mutsig side bar plot
# basal_comut.axes['Mutation_frequency'].axvline(4, color = 'black', linestyle = 'dotted', linewidth = 1)

### will get back to this later
basal_comut.axes['Mutation_frequency'].axvline(7, color = 'black', linestyle = 'dotted', linewidth = 1)
basal_comut.axes['Mutation_frequency'].axvline(13, color = 'black', linestyle = 'dotted', linewidth = 1)
basal_comut.axes['Mutation_frequency'].set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
                                                   10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
                                                   20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5,
                                                   30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5,
                                                   40.5, 41.5, 42.5, 43.5])
basal_comut.axes['Mutation_frequency'].set_yticklabels(list(percentages))
basal_comut.axes['Mutation_frequency'].tick_params(axis='y', pad=-34)



# basal_comut.axes['CN_frequency'].axvline(6, color = 'black', linestyle = 'dotted', linewidth = 1)
# basal_comut.axes['CN_frequency'].axvline(12, color = 'black', linestyle = 'dotted', linewidth = .75)

mut_type_leg = basal_comut.add_axis_legend(name='Event Type', bbox_to_anchor = (1, 0.89), title='Event Type', order = mut_order)
# basal_comut.figure.savefig('/fh/fast/ha_g/projects/NelsonLab/MedGenome_WGS_LNCaP/robinson_comut_figures/tpersseComutWork/plots/lncapComut_v1_20241231_allSamples.NoSRR.pdf', format='pdf', bbox_inches = 'tight', dpi = 300)

basal_comut.figure.savefig('outputComut.pdf', format='pdf', bbox_inches = 'tight', dpi = 300)
