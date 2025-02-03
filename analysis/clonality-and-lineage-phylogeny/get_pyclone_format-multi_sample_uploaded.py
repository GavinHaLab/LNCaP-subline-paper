#Import the necessary libraries
import numpy as np
import pandas as pd
import io
import os
import vcf
import pathlib
import glob
from pathlib import Path
import sys
import matplotlib.pyplot as plt

__Purpose__ = "#SCRIPT to use the generated GATK files to create the PyClone input files based on the different copy number inputs from TitanCNA"

try:
    import pandasql as ps
    # print("module 'mutagen' is installed")
except ModuleNotFoundError:
    print("module 'pandasql' is not installed")
    # or
    pip.main(['install', pandasql]) # the install function from the question
# import pandasql as ps

def read_vcf(path, sample):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str, sample: str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

## SET UP DIRECTORIES
curdir = Path(os.getcwd())
project_path = str(curdir.parents[2]) + '/'

print("Current Main Directory: ", project_path)
print("")


mutation_dir = "/pathtodir/mutation_calling/active/mutect2_06162023/" # specifiy the path to mutect2 or your mutation calling output files
titan_dir = "/pathtodir/TitanCNA_repaired/results/" # specify the path to TitanCNA or your copy number caller output files
gatk_dir = "/pathtodir/pyclone/multi_sample/GATK_results/" # specify the path to the GATK files or the directory containing read counts

output_dir = "/pathtodir/mpyclone/multi_sample/inputs/"
print("Current Output Directory: ", output_dir)

# Get naming conventions to correlate mutect2 with titanCNA
all_mut_entries = sorted(os.listdir(str(mutation_dir)))
all_cna_entries = sorted(os.listdir(str(titan_dir)))
mutation_directories = ([entry for entry in all_mut_entries if os.path.isdir(os.path.join(mutation_dir, entry))])
cna_directories = ([entry for entry in all_cna_entries if os.path.isdir(os.path.join(titan_dir, entry))])
samples_dict = dict(zip(mutation_directories, cna_directories))
# print("Sample Naming Conventions:   ", samples_dict)

print("     ...Finished reading files and setting up directories")

def alt_count_splits(row, sample):
    alt_count = row[sample].split(':')[1].split(',')[1]
    return alt_count

def ref_count_splits(row, sample):
    ref_count = row[sample].split(':')[1].split(',')[0]
    return ref_count


desired_samples = [['LNCaP_FGC', 'LNCaP_clone_FGC', 'LNCaP_16D', 'LNCaP_42D', 'LNCaP_42F', 'LNCaP_C4', 'LNCaP_C42', 'LNCaP_C42B', 'LNCaP_Abl', 'LNCaP_95', 'LNCaP_APIPC', 'LNCaP_shAR', 'LNCaP-AR-907', 'LNCaP-AR-909','LNCaP_FGC_SRR7943697']]

missing = []
count = 1
## Loop through all of the samples
for patient in desired_samples:
    for s in patient:   
        print("Processing sample: ", s)
        # Load in GATK data per sample
        gatk_file = pd.read_csv(gatk_dir + s + '_allelicCounts.tsv')
        gatk_file_contig = gatk_file.loc[(gatk_file['@HD\tVN:1.6']).str.contains('CONTIG')]
        gatk_file2 = (gatk_file.iloc[(gatk_file_contig.index.tolist()[0]):]).reset_index(drop=True)
        gatk_file2.columns = gatk_file2.iloc[0]
        gatk_file2.drop(gatk_file2.index[0], inplace=True)        
        mut_vcf1 = gatk_file2['CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE'].str.split('\t', expand=True)
        mut_vcf1.columns = ['CONTIG', 'POSITION', 'REF_COUNT', 'ALT_COUNT', 'REF_NUCLEOTIDE', 'ALT_NUCLEOTIDE']

        # Filter out indels
        mut_vcf = (mut_vcf1.loc[(mut_vcf1['REF_NUCLEOTIDE'].str.len()) == mut_vcf1['ALT_NUCLEOTIDE'].str.len()]).reset_index(drop=True)

        # Add to new dataframe
        py_vcf = pd.DataFrame()
        py_vcf['mutation_id'] = mut_vcf['CONTIG'] + ':' + mut_vcf['POSITION'].astype(str) + '-' + mut_vcf['POSITION'].astype(str)
        py_vcf['chr'] = mut_vcf['CONTIG']
        py_vcf['start'] = mut_vcf['POSITION']
        py_vcf['end'] = mut_vcf['POSITION']
        py_vcf['REF'] = mut_vcf['REF_NUCLEOTIDE']
        py_vcf['ALT'] = mut_vcf['ALT_NUCLEOTIDE']

        # Get ref and alt counts
        # py_vcf['ref_counts']= mut_vcf.apply(lambda x: ref_count_splits(x, s), axis=1)
        # py_vcf['alt_counts']= mut_vcf.apply(lambda x: alt_count_splits(x, s), axis=1)
        py_vcf['ref_counts'] = mut_vcf['REF_COUNT']
        py_vcf['alt_counts'] = mut_vcf['ALT_COUNT']

        # drop indels
        py_vcf = py_vcf.drop(py_vcf[py_vcf['REF'] == '-'].index)
        py_vcf = py_vcf.drop(py_vcf[py_vcf['ALT'] == '-'].index)
        py_vcf.reset_index(drop=True)

        # drop the REF and ALT columns
        py_vcf.drop(columns=['REF', 'ALT'], inplace=True)

        cnv_file = glob.glob(''.join([titan_dir, s, '*.titan.ichor.seg.noSNPs.txt']))
        somatic_cnv = pd.read_csv(cnv_file[0], sep='\t')
        print("     ...Finished reading CNV data for sample ", s)

        sqlcode = ''' 
        select py_vcf.*, somatic_cnv.Corrected_Copy_Number, somatic_cnv.Corrected_MajorCN, somatic_cnv.Corrected_MinorCN
        from py_vcf
        left join somatic_cnv
        on py_vcf.chr = somatic_cnv.Chromosome
        where py_vcf.start >= somatic_cnv.Start and py_vcf.end <= somatic_cnv.End
        '''

        final_df = ps.sqldf(sqlcode,locals())

        # Mostly to check for chrX values not working properly
        gen = (x for x in list(final_df.index.values) if ((final_df.iloc[x])['chr'] == 'chrX'))
        for x in gen:
            if (pd.isna(final_df.at[x, 'Corrected_MajorCN'])):
                final_df.at[x, 'Corrected_MajorCN'] = final_df.at[x, 'Corrected_Copy_Number']
                final_df.at[x, 'Corrected_MinorCN'] = 0

        final_df.drop(['chr', 'start', 'end'], axis=1, inplace=True)

        # print(final_df)

        ####final_df = final_df.rename(columns={'Corrected_Copy_Number': 'normal_cn', 'Corrected_MajorCN': 'major_cn', 'Corrected_MinorCN': 'minor_cn'}) # Incorrect line normal; CN should be 2
        final_df = final_df.rename(columns={ 'Corrected_MajorCN': 'major_cn', 'Corrected_MinorCN': 'minor_cn'})
        final_df['normal_cn'] = 2
        # drop correct copy number columns
        
        # fill in all blanks with NaN
        final_df.fillna(0, inplace=True)
        print("     ...Getting purity information")
        ## Get tumour_content
        purity = []
        params_file = glob.glob(''.join([titan_dir, s, '*.params.txt']))
        print("     ...Reading purity data for sample ", params_file)
        # file = titan_dir + '/' + '{0}.params.txt'.format(cn_sample)
        text_file = open(params_file[0], "r")

        #read whole file to a string
        data = text_file.read()
        #close file
        text_file.close()
        strs = data.split('\t')
        # Purity
        purity_data = (1 - pd.to_numeric(strs[1].split('\n')[0]))
        tumour_content_series = [purity_data] * len(final_df['mutation_id'])
        final_df['tumour_content'] = tumour_content_series

        s_series = ['{0}'.format(s)] * len(final_df['mutation_id'])
        final_df.insert(1, 'sample_id', s_series)

        ## Round copy numbers
        final_df['normal_cn'] = (final_df['normal_cn'].round())
        final_df['major_cn'] = (final_df['major_cn'].round())
        final_df['minor_cn'] = (final_df['minor_cn'].round())
        
        # final_df.fillna(0, inplace=True)
        ##final_df = final_df.astype({"normal_cn": "int", "major_cn": "int", "minor_cn": "int"}) #unsure why this is here
        # plot VAF vs minor_cn
        # esnure all columns are numeric
        final_df['ref_counts'] = pd.to_numeric(final_df['ref_counts'])
        final_df['alt_counts'] = pd.to_numeric(final_df['alt_counts'])

        final_df['VAF'] = final_df['alt_counts'] / (final_df['ref_counts'] + final_df['alt_counts'])
        final_df.plot.scatter(y='VAF', x='minor_cn')
        out_file = output_dir + s + '_VAf_minorCN.png'
        plt.savefig(out_file)
        # drop the column VAF
        final_df.drop(columns="VAF", inplace=True)
        # final_df.to_csv(output_dir + s + '.tsv', sep='\t', index=False)
        print("     ...Completed PyClone format generation for sample: ", s)
        print('')
        count += 1
    
        # Continue aggregating if it is not the first sample in patient group
        if (s == patient[0]):
            output_df = final_df.copy()
        else:
            output_df = pd.concat([output_df, final_df], ignore_index=True)

    ##output_df.to_csv(output_dir + 'Subline_' + patient[0] + '.tsv', sep='\t', index=False)
    output_df.to_csv(output_dir + 'All14lines_LNCaP.tsv', sep='\t', index=False)
    print("\t...Completed PyClone format generation for patient: ", s)
    
print("Completed for all samples!")

# retain only the same mutation IDs across all samples
num_samples = len(desired_samples[0])
output_df = output_df[output_df['mutation_id'].isin(output_df['mutation_id'].value_counts()[output_df['mutation_id'].value_counts() == num_samples].index)]
# write to file
output_df.to_csv(output_dir + 'All14lines_LNCaP_commonMutations.tsv', sep='\t', index=False)
print("the final length of the dataframe with unique mutations is: ", len(output_df))

# post processing needed to fix the copy number values and column names
#drop the column corrected_copy_number
output_df.drop(columns="Corrected_Copy_Number", inplace=True)
# If mutation id strats with chrX, then set normal_cn to 1
output_df.loc[output_df['mutation_id'].str.startswith('chrX'), 'normal_cn'] = 1
# MAKE SURE THESE ARE INTEGERS
output_df['normal_cn'] = output_df['normal_cn'].astype(int)
output_df['major_cn'] = output_df['major_cn'].astype(int)
output_df['minor_cn'] = output_df['minor_cn'].astype(int)
output_df.to_csv(output_dir + 'All14lines_LNCaP_commonMutations_fixed.tsv', sep='\t', index=False)
