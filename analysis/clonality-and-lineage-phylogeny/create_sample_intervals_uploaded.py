# Scripts to generate intervals for each sample and then combine them to create a final interval file to obtain read counts for WT and ALT allele for each sample from the bam file using GATK CollectAllelicCounts tool.


import numpy as np
import pandas as pd
import io
import os
import sys
import matplotlib.pyplot as plt


mut_dir = '/pathtodir//tumor_only/results/' # path to mutect2 or mutation calling output files


# list of samples   
samples = [['LNCaP_FGC', 'LNCaP_clone_FGC', 'LNCaP_16D', 'LNCaP_42D', 'LNCaP_42F', 'LNCaP_C4', 'LNCaP_C42', 'LNCaP_C42B', 'LNCaP_Abl', 'LNCaP_95', 'LNCaP_APIPC', 'LNCaP_shAR', 'LNCaP-AR-907', 'LNCaP-AR-909']]
for patient in samples:
    patient_intervals = []
    for sample in patient:
        print('Processing sample: ', sample)
        ## Get sample data - mutations for UC
        muts = pd.read_csv(mut_dir + '/' + sample + '/pass_snvs.annovar.hg38_multianno.txt', sep='\t')
        # added  to exlcude benign variants from clinvar and retain potentially pathogenic variants
        functional_df = muts.loc[
        (muts['cosmic70'] != '.') |
        (muts['SIFT_pred'] == 'D') |
        (muts['LRT_pred'] == 'D') |
        (muts['MutationTaster_pred'] == 'D') |
        (muts['FATHMM_pred'] == 'D') |
        (~muts['CLNSIG'].isin(['.', 'Likely_benign', 'Benign', 'Benign/Likely_benign', 'not_provided']))
    ]
        ##functional_df = muts # uncomment this line to get all mutations
        ############################################# POTENTIAL CHANGES END #############################################

        ### Remove Common Variants - GNOMAD ####
        ##threshold = 0.1
        threshold = 0.01 #this is 1% threshold
        functional_df.loc[functional_df['gnomAD_genome_ALL'] == '.', 'gnomAD_genome_ALL'] = '0'
        out = functional_df[pd.to_numeric(functional_df['gnomAD_genome_ALL']) > threshold].index
        filtered_df_init = functional_df.drop(out)
        #### Remove Common Variants - EXAC ####
        filtered_df_init.loc[filtered_df_init['ExAC_ALL'] == '.', 'ExAC_ALL'] = '0'
        out = filtered_df_init[pd.to_numeric(filtered_df_init['ExAC_ALL']) > threshold].index
        filtered_df = filtered_df_init.drop(out)

        #### GET VAF ####
        # For N-T Mode
        info_fields = filtered_df[['Otherinfo12', 'Otherinfo13']]
        vafs = []
        vaf_index = list(filtered_df.index)
        ads_ref = []
        ads_alt = []
        dps = []
        for x in vaf_index:
            string_fields = info_fields['Otherinfo13'][x]
            sf = string_fields.split(':')
            a_d_ref = sf[1].split(',')[0]
            a_d_alt = sf[1].split(',')[1]
            dp = sf[3]
            vaf = pd.to_numeric(a_d_alt) / pd.to_numeric(dp)
            ads_ref.append(a_d_ref)
            ads_alt.append(a_d_alt)
            dps.append(dp)
            vafs.append(vaf)
        vaf_col = pd.Series(vafs, name='VAF').set_axis(vaf_index, inplace=False)
        ad_alt_col = pd.Series(ads_alt, name='AD_ALT').set_axis(vaf_index, inplace=False)
        ad_ref_col = pd.Series(ads_ref, name='AD_REF').set_axis(vaf_index, inplace=False)
        dp_col = pd.Series(dps, name='DP').set_axis(vaf_index, inplace=False)
        filtered_df.insert(13, 'VAF', vaf_col)
        filtered_df.insert(14, 'AD_ALT', ad_alt_col)
        filtered_df.insert(15, 'AD_REF', ad_ref_col)
        filtered_df.insert(16, 'DP', dp_col)

        #### FINAL: Filter by VAF/Read Counts ######
        # trial 1 was the original approach and was later modified to the current approach i.e.trial2 which is more stringent
        # Filter by VAF
        f_vaf_Trial1 = filtered_df.loc[filtered_df['VAF'] >= 0.1] # original trial 1
        # retain all mutations whose VAf is between 0 - 0.45 and between 0.60 to 0.80
        f_vaf = filtered_df.loc[((filtered_df['VAF'] >= 0.1) & (filtered_df['VAF'] <= 0.45)) | ((filtered_df['VAF'] >= 0.60) & (filtered_df['VAF'] <= 0.80))] # trial 2 cut offs
        # rescue those mutations from filtered that do not fit the above criteria but are in COSMIC
        f_vaf = pd.concat([f_vaf, filtered_df.loc[(filtered_df['cosmic70'] != '.') & ((filtered_df['VAF'] >= 0.80) | ((filtered_df['VAF'] > 0.45) & (filtered_df['VAF'] < 0.60)))]]) # trial2 cut offs

        ## Filter by total count in tumor
        f_DP = f_vaf.loc[pd.to_numeric(f_vaf['DP']) >= 10] # trial 2 cut offs
        f_DP_Trial1 = f_vaf_Trial1.loc[pd.to_numeric(f_vaf_Trial1['DP']) >= 10] # original trial 1
        ## Filter by total alternate reads
        f_ALT_COUNT = f_DP.loc[pd.to_numeric(f_DP['AD_ALT']) >= 3] # # trial 2 cut offs
        f_ALT_COUNT_Trial1 = f_DP_Trial1.loc[pd.to_numeric(f_DP_Trial1['AD_ALT']) >= 3]
        ## Filter by total normal reads
        # f_REF_COUNT = f_ALT_COUNT.loc[pd.to_numeric(f_ALT_COUNT['AD_REF']) >= 3]
        f_final = (f_ALT_COUNT.copy()).reset_index(drop=True) # trial 2 cut offs
        f_final_Trial1 = (f_ALT_COUNT_Trial1.copy()).reset_index(drop=True) # original trial 1

        print('     ............Size of Final Mutations for sample ', sample, ':  ', f_final.shape)
        sample_interval = (f_final['Chr'] + ':' + f_final['Start'].astype(str) + '-' + f_final['End'].astype(str)).tolist()
        patient_intervals = patient_intervals + sample_interval
        print("     ............Size of current patient intervals: ", len(patient_intervals))
        # write each sample to a file
        out_file = mut_dir + '/' + sample + '/' + sample + '_functionalMutations_Trial2.txt'
        f_final.to_csv(out_file, sep='\t', index=False)
        # plot a histogram of VAF values
        f_final['VAF'].plot.hist(bins=20, alpha=0.5, title=sample,edgecolor='black')
        # axis labels
        out_vaf_plot = mut_dir + '/' + sample + '/' + sample + '_VAF_hist.png'
        plt.savefig(out_vaf_plot)
        # clear the plot
        plt.clf()
 
        # write the originam Trial1s approach to a file
        out_file_Trial1 = mut_dir + '/' + sample + '/' + sample + '_functionalMutations_Trial1.txt'
        f_final_Trial1.to_csv(out_file_Trial1, sep='\t', index=False)
        # plot a histogram of VAF values
        f_final_Trial1['VAF'].plot.hist(bins=20, alpha=0.5, title=sample,edgecolor='black')
        out_vaf_plot_Trial1 = mut_dir + '/' + sample + '/' + sample + '_VAF_hist_Trial1.png'
        plt.savefig(out_vaf_plot_Trial1)
        plt.clf()
    
    # After all samples per group are finished, remove duplicates and sort
    patient_intervals = list(set(patient_intervals))
    patient_intervals.sort()

    #output = os.getcwd() + '/LNCaP_all_updated1Sample_deleterious.intervals'
    output = os.getcwd() + '/LNCaP_14samples_functionalMutations_secondSet.intervals'
    print("     ............Final size of patient intervals: ", len(patient_intervals))
    with open(output, "w") as file:
        # Write each element of the list to the file
        for item in patient_intervals:
            file.write(str(item) + "\n")

    print("")