#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import multiprocessing as mp


#### Directories ####
samp_dir = "/path/to/mutect/results/" ## directory containing all the samples
# indel_dir = "/fh/fast/ha_g/projects/NelsonLab/MedGenome_WGS_LNCaP/TitanCNA/TitanCNA_SV_WGS_tumorOnly_10kb/results/svaba_tumor_only/"
out_dir = "/path/to/output/directory/" ## output directory
maneFile = "path/to/MANE.GRCh38.v1.4.summary.txt" ## get path for mane summary, use this to get transcript ID...

## load in mane summary file, convert to dict with key = gene, subkey = mane_select, mane_plus_clinical, values = transcript ID

mane_df = pd.read_csv(maneFile, sep='\t', header=0) ## ensure it's tab separated
mane_dict = {}
# set gene col as index
mane_df.set_index('symbol', inplace=True)
for gene in mane_df.index:
    mane_dict[gene] = {}
    ## source column for mane_select and mane_plus_clinical is MANE_status. Note that there are spaces in the values( 'MANE Select' and 'MANE Plus Clinical') so will want to replace these with underscores
    # check if multiple transcripts are present for a gene
    if isinstance(mane_df.loc[gene, 'RefSeq_nuc'], str):
        if mane_df.loc[gene, 'MANE_status'] == 'MANE Select':
            mane_dict[gene]['mane_select'] = mane_df.loc[gene, 'RefSeq_nuc']
        elif mane_df.loc[gene, 'MANE_status'] == 'MANE Plus Clinical':
            mane_dict[gene]['mane_plus_clinical'] = mane_df.loc[gene, 'RefSeq_nuc']

    else:
        if 'MANE Select' in mane_df.loc[gene, 'MANE_status'].values:
            select_transcript = mane_df.loc[gene, 'RefSeq_nuc'][mane_df.loc[gene, 'MANE_status'] == 'MANE Select']
            mane_dict[gene]['mane_select'] = select_transcript.values[0]
        elif 'MANE Plus Clinical' in mane_df.loc[gene, 'MANE_status'].values:
            clinical_transcript = mane_df.loc[gene, 'RefSeq_nuc'][mane_df.loc[gene, 'MANE_status'] == 'MANE Plus Clinical']
            mane_dict[gene]['mane_plus_clinical'] = clinical_transcript.values[0]

selectCt = 0
clinicalCt = 0
for key in mane_dict:
    print(key, mane_dict[key])
    if 'mane_select' in mane_dict[key]:
        selectCt += 1
    if 'mane_plus_clinical' in mane_dict[key]:
        clinicalCt += 1

# print(selectCt, clinicalCt)
## get count of subkeys 


samp_full_list = os.listdir(samp_dir)
# remove samples with clone in the name
# samp_full_list = [s for s in samp_full_list if 'clone' not in s]
print(samp_full_list)

def processFile(s):
    print(s)
    # if s == 'LNCaP_shAR'
    #### Get sample - original dataframe ####
    file = samp_dir+s+'/pass_variants.annovar.hg38_multianno.txt' ## changed from pass_variants... we only want snvs from mutect results
    # indel = indel_dir+s+f'/{s}.svaba.indel.hg38_multianno.txt'
    if s == 'LNCaP_shAR':
        s = s+'-pATK' ## some naming convention changes here...

    if os.path.exists(file):
        print('file exists')
        # print(file, indel)
        df = pd.read_csv(file, sep='\t')
        # df_indel = pd.read_csv(indel, sep='\t')

        # for easier handling, will concatenate the two dataframes after adding 'mutect' and 'svaba' as 'type' column
        df['Source'] = 'mutect'

        ## add parsing here to figure out if variant is snv or indel
        ## look at ref and alt. if len(ref) == len(alt) == 1 and they are in ['A', 'T', 'C', 'G'], then it's an snv, otherwise it's an indel
        df['Type'] = 'SNV'
        df.loc[(df['Ref'].str.len() > 1) | (df['Alt'].str.len() > 1), 'Type'] = 'INDEL' ## if either ref or alt is greater than 1, then it's an indel

        ## should also check if both ref and alt are in ['A', 'T', 'C', 'G'] for those that are still SNVs
        df.loc[~df['Ref'].isin(['A', 'T', 'C', 'G']) | ~df['Alt'].isin(['A', 'T', 'C', 'G']), 'Type'] = 'INDEL'

        # df_indel['Source'] = 'svaba'
        # df = pd.concat([df, df_indel])
        # print(df.head)
        # print([x for x in df.columns])
        ## need to reindex the dataframe
        # df.reset_index(drop=True, inplace=True)
        ## need to add
        ## should concatenate if possible...
        #### Remove Common Variants ####
        threshold = 0.1
        cols=df.columns
        # df.loc[df['AF'] == '.', 'AF'] = '0' # here, are we trying to remove gnomad variants?...
        df.loc[df['gnomAD_genome_ALL'] == '.', 'gnomAD_genome_ALL'] = '0'
        df.loc[df['ExAC_ALL'] == '.', 'ExAC_ALL'] = '0'

        # now, filter out variants with AF > threshold
        
        out = df[pd.to_numeric(df['gnomAD_genome_ALL']) > threshold].index
        filtered_df = df.drop(out)
        out = df[pd.to_numeric(df['ExAC_ALL']) > threshold].index
        filtered_df = df.drop(out)


        #### Get Cosmic Variants ####
        cosmics = filtered_df.loc[filtered_df['cosmic70'] != '.']
        cosmic_column = []
        # print(cosmics['cosmic70'])

        cosmic_index = list(filtered_df.index) 

        for row in cosmic_index:
            # print(s, row)
            # here, we are trying to get the cosmic70 ID for each variant... will use this to filter out variants that are not in cosmic
            val = filtered_df.loc[row]
            geneval = val['Gene.refGene']
            # print(geneval)
            # print(cosmics['Gene.refGene'].values)
            if geneval in cosmics['Gene.refGene'].values: # here, we are checking if the gene is in the cosmic70 dataframe... 
                cosID_df = filtered_df.loc[row]['cosmic70']
                cosmic_column.append(cosID_df)
            else:
                cosmic_column.append('.')
        # print(cosmic_column)
        ## here, we are adding the cosmic70 column to the dataframe
        cos_col = pd.Series(cosmic_column, name='Cosmic70').set_axis(cosmic_index).rename('Cosmic70') # create series with the cosmic70 column, and set the index to the cosmic_index (which is the index of the filtered dataframe)
        ## now, we insert the cosmic70 column into the dataframe
        filtered_df.insert(9, 'Cosmic70', cos_col) 

        #### Get Intervar variants ####
        intervars = filtered_df.loc[filtered_df['InterVar_automated'] != '.']
        intervar_column = []
        int_index = list(filtered_df.index)
        for row in int_index:
            val = filtered_df.loc[row]
            geneval = val['Gene.refGene']
            if geneval in intervars['Gene.refGene'].values:
                intervar_ID = filtered_df.loc[row]['InterVar_automated']
                intervar_column.append(intervar_ID)
            else:
                intervar_column.append('.')

        int_col = pd.Series(intervar_column, name='InterVar').set_axis(int_index).rename('InterVar')
        filtered_df.insert(10, 'InterVar', int_col)

        #### Get Clinvar variants ####
        clinvars = filtered_df.loc[filtered_df['CLNSIG'] != '.']
        clinvar_column = []
        cln_index = list(filtered_df.index)

        for row in cln_index:
            val = filtered_df.loc[row]
            geneval = val['Gene.refGene']
            if geneval in clinvars['Gene.refGene'].values:
                clinvar_ID = filtered_df.loc[row]['CLNSIG']
                clinvar_column.append(clinvar_ID)
            else:
                clinvar_column.append('.')

        clin_col = pd.Series(clinvar_column, name='ClinVar').set_axis(cln_index).rename('ClinVar')
        filtered_df.insert(11, 'ClinVar', clin_col)

        #### GET VAF ####
        info_fields = filtered_df[['Otherinfo12', 'Otherinfo13']]
        vafs = []
        vaf_index = list(filtered_df.index)
        ads_ref = []
        ads_alt = []
        dps = []

        for x in vaf_index:
            string_fields = info_fields['Otherinfo13'][x]
            sf = string_fields.split(':')
            if filtered_df.loc[x, 'Source'] == 'mutect':
                a_d_ref = sf[1].split(',')[0]
                a_d_alt = sf[1].split(',')[1]
                dp = sf[3]
                vaf = pd.to_numeric(a_d_alt) / pd.to_numeric(dp)
                ads_ref.append(a_d_ref)
                ads_alt.append(a_d_alt)
                dps.append(dp)
                vafs.append(vaf)
            else:
                a_d_alt = sf[1]
                dp = sf[2]
                vaf = pd.to_numeric(a_d_alt) / pd.to_numeric(dp)
                a_d_ref = pd.to_numeric(dp) - pd.to_numeric(a_d_alt)
                ads_ref.append(a_d_ref)
                ads_alt.append(a_d_alt)
                dps.append(dp)
                vafs.append(vaf)
                print(a_d_ref, a_d_alt, dp, vaf) ## check if this is working correctly
        


        vaf_col = pd.Series(vafs, name='VAF').set_axis(vaf_index)
        ad_alt_col = pd.Series(ads_alt, name='AD_ALT').set_axis(vaf_index)
        ad_ref_col = pd.Series(ads_ref, name='AD_REF').set_axis(vaf_index)
        dp_col = pd.Series(dps, name='DP').set_axis(vaf_index)
        filtered_df.insert(13, 'VAF', vaf_col)
        filtered_df.insert(14, 'AD_ALT', ad_alt_col)
        filtered_df.insert(15, 'AD_REF', ad_ref_col)
        filtered_df.insert(16, 'DP', dp_col)

        # set cols to be int values
        filtered_df['DP'] = pd.to_numeric(filtered_df['DP'])
        filtered_df['AD_ALT'] = pd.to_numeric(filtered_df['AD_ALT'])
        filtered_df['AD_REF'] = pd.to_numeric(filtered_df['AD_REF'])
        filtered_df['VAF'] = pd.to_numeric(filtered_df['VAF'])
        # filter out variants with VAF < 0.1
        filtered_df = filtered_df[filtered_df['VAF'] > 0.1]
        # filter out variants with DP < 10 and AD_ALT < 3
        filtered_df = filtered_df[filtered_df['DP'] > 10]
        filtered_df = filtered_df[filtered_df['AD_ALT'] > 3]

        #### Create final DF with changed column names ####
        final_df = filtered_df[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 
                                'ExonicFunc.refGene', 'Cosmic70', 'InterVar', 'ClinVar', 'cytoBand', 'AD_REF', 'AD_ALT', 'DP', 'VAF', 'AAChange.refGene', 'Source', 'Type']]
        #### Re-format column names for easier parsing ####
        final_interest_df = final_df.rename(columns={'Chr': 'CHR', 'Ref': 'REF', 'Alt': 'ALT', 'Gene.refGene': 'Gene',
                                                        'ExonicFunc.refGene': 'Mutation_Type', 'AAChange.refGene': 'AA_Change'})

        #### Convert sample name ####
        sf = [i for i in range(len(samp_full_list)) if samp_full_list[i] == s] # get the index of the sample in the sample list
        sample = s
        samp_id = ['{0}'.format(sample)] * len(final_interest_df['Gene']) # create a list of the sample name
        final_interest_df.insert(0, 'Sample', samp_id) # insert the sample name into the dataframe

        #### Splice the FIRST AA_Change ####
        aa_index = list(final_interest_df.index)
        all_aa = []

        ## we want to use the amino acid change which corresponds to the mane_select transcript if it exists, otherwise use the mane_plus_clinical transcript if it exists, otherwise use the first transcript in the list
        # easiest way to do this is to take values from AAchange.refGene, split on comma, store element in dict where key is transcript ID, value is other info joined by commas. If any key found in mane_dict, use that value, otherwise use the first value in the dict
        aa_changes = final_interest_df['AA_Change']
        # aa_dict = {}
        # need to get length of aa_changes, figure out
        for i in aa_index:
            aa_dict = {}
            aa_change = aa_changes[i]
            if aa_change == '.':
                # aa_dict[i] = '.'
                all_aa.append('.')
            elif aa_change == 'UNKNOWN':
                # aa_dict[i] = 'UNKNOWN'
                all_aa.append('UNKNOWN')
            else:
                aa_change = aa_change.split(',')
                for change in aa_change:
                    transcript = change.split(':')[0]
                    if transcript not in aa_dict:
                        aa_dict[transcript] = ":".join(change.split(':')[2:5])
                ## now check if mane_select or mane_plus_clinical is in the dict
                gene = final_interest_df.loc[i, 'Gene']
                if gene in mane_dict:
                    if mane_dict[gene]['mane_select'] in aa_dict: # if mane_select transcript is in the dict
                        # return the value of of this key from aa_dict
                        all_aa.append(aa_dict[mane_dict[gene]['mane_select']])
                    # elif mane_dict[gene]['mane_plus_clinical'] in aa_dict:
                    elif 'mane_plus_clinical' in mane_dict[gene]:
                        if mane_dict[gene]['mane_plus_clinical'] in aa_dict:
                            all_aa.append(aa_dict[mane_dict[gene]['mane_plus_clinical']])
                    else:
                        ## switch to using the first value in the dict
                        all_aa.append(list(aa_dict.values())[0])
                else:
                    # all_aa.append(aa_dict[sorted(aa_dict.keys())[0]]) # throwing error here: TypeError: '<' not supported between instances of 'str' and 'int'
                    all_aa.append(list(aa_dict.values())[0])
        f = final_interest_df.drop(columns='AA_Change')
        f['AA_Change'] = all_aa

        return f ## think about what to return here
    else:
        print(f'file for {s} does not exist')
        # return pd.DataFrame()
## create a pool of workers
pool = mp.Pool(15)
results = pool.map(processFile, samp_full_list)
pool.close()
pool.join()

all_samples_df = pd.concat(results)

all_samples_df.to_csv(out_dir+'/LNCaP_WGS_TumorOnly_allVariants_Manifest_20241213_mutectOnly.txt', index=False, sep='\t')



