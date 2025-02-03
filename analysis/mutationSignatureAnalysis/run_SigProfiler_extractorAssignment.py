#!/usr/bin/python
__author__="Pushpa Itagi"
__purpose__="To generate a SNV matrix and then fit this matrix to COSMIC SBS signatures"

print("loading modules")
import SigProfilerMatrixGenerator
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerMatrixGenerator.scripts import SVMatrixGenerator as sv
import pandas as pd

import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
import sys
import os
print(sys.path)

def main_function():
    print(" In sig profiler matrix generator Function")
    samples = ['LNCaP_all_excludes']
    for s in samples:
        input_file_dir = "/path/to/sigProfiler/input"

        project_name = "{0}".format(s)
        matrices = matGen.SigProfilerMatrixGeneratorFunc(project_name, "GRCh38",input_file_dir, plot=True, exome=False, bed_file=None,tsb_stat=False, seqInfo=True, cushion=100)
    return 0;

def sigProfilerAssignment():
    #set directories and paths to signatures and samples
    print(" In sig profiler sigProfilerAssignment  Function")

    dir_inp = "/path/to/sigProfiler/input"

    samples_file     = dir_inp+"/output/SBS/LNCaP_all_excludes.SBS96.all"      # FOR SBS
    output_dir      =  "/path/to/sigProfiler/output/SBS" #Output Directory
    sigs        = "COSMIC_v3_SBS_noSBS42_GRCh38.txt" #Custom Signature Database, IF NEEDED ^^^^^

    exclude_groups = ['UV_signatures', 'Artifact_signatures']

    # Analysis of SP Assignment
    Analyze.cosmic_fit(samples_file,
                    output_dir,
                    genome_build="GRCh38",
                    cosmic_version=3,
                    verbose=True,
                    context_type="96",
                    input_type="matrix",
                    collapse_to_SBS96=True,
                    make_plots=True,
                    signature_database=sigs,
                    exclude_signature_subgroups=exclude_groups,
                    sample_reconstruction_plots=True,
                    exome=False,export_probabilities=True)


if __name__=="__main__":
    main_function()
    sigProfilerAssignment()

