#!/usr/bin/env python3

import os, argparse

def get_args():
    parser = argparse.ArgumentParser(description="Obtain purity and ploidy from TitanCNA output")
    parser.add_argument("-i", "--indir", help="directory upstream of titan results dirs", required=True)
    # parser.add_argument("-o", "--outdir", help="files t", required=True)
    return parser.parse_args()

args = get_args()

if not args.indir.endswith('/'):
    args.indir += '/'


outfile = 'purity_ploidy.txt'


with open(outfile, 'w') as o:
    o.write('Sample\tPurity\tPloidy\n')
    for file in os.listdir(args.indir):
        if file.endswith('.params.txt'):
            with open(args.indir + '/' + file, 'r') as f:
                for line in f:
                    if line.startswith('Normal contamination estimate:'):
                        purity = str(1 - float(line.split('\t')[1].strip()))
                    if line.startswith('Average tumour ploidy estimate:'):
                        ploidy = line.split('\t')[1].strip()
            name = file.split('.params.txt')[0]
            o.write(name + '\t' + purity + '\t' + ploidy + '\n')




    