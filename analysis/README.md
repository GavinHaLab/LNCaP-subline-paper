# Code and scripts for analyses associated with the manuscript.

Scripts used for each analysis were grouped by analysis type

## comut

<strong>generateComut.py</strong>

Script used to read in mutation data formatted by makeComutInput.py, manipulate into proper format and combine with CN calls, then produce comut plots seen in figures 1E, 2A, S2B

<strong>getMutFreq.py</strong>

Script used to pull mutation frequency information from input comut-formatted SNV/Indel calls, for use in comut figures 

<strong>makeComutInput.py</strong>

Script used to pull relevant SNV/Indel calls from manifest file, reformat to proper comut input format and subset to genes of interest

<strong>reformatCNdata.py</strong>

Script to do slight adjustments to outputs from generate_comut_CN.R (see makeMatrixFromTitan section) to match sample naming schemes


## lohScore

<strong>cal_loh_score_cp_adjusted.R</strong>

Script used to calculate LOH scores seen in table 3 from TitanCNA outputs

<strong>run.sh</strong>

Example bash script showing how to run cal_loh_score_cp_adjusted.R

## makeMatrixFromTitan

<strong>generate_comut_CN.R</strong>

Rscript used to reformat makeMatrixFromTITAN-ICHOR_segBased_hg38.R output segment-based gene-level CN matrix to comut input format, as well as providing several thresholded values for CN (thresholds based on CP)

<strong>makeMatrixFromTITAN-ICHOR_segBased_hg38.R</strong>

Rscript to pull gene-level copy number, LOH, etc. calls from TitanCNA outputs and produce several matrices with genes as columns and samples as rows

<strong>make_samplelist.py</strong>

Script to produce a mandatory input for makeMatrixFromTITAN, purity_ploidy.txt, given a Titan run directory. Output file is tab-delimited file containing sample name, estimated purity value, and estimated ploidy value

<strong>runTitanMatrix.sh</strong>

Example bash script showing how to run makeMatrixFromTITAN-ICHOR_segBased_hg38.R

## mutationSignatureAnalysis

<strong>COSMIC_v3_SBS_noSBS42_GRCh38.txt</strong>

Input SBS matrix with SBS42 removed, used as input when running run_SigProfiler_extractorAssignment.py. Outputs were mutation signature matrices, which were used to produce values seen in table 3, as well as plots seen in figures S1F, S2H, and S2I

<strong>makeInputVcfs.py</strong>

Script used to format and filter input vcfs directly from mutect2 run dir, to be used in run_SigProfiler_extractorAssignment.py 

<strong>run_SigProfiler_extractorAssignment.py</strong>

Script used to run sigProfilerMatrixGenerator and SigProfilerAssignment on all LNCaP subline samples

## snvAndIndelManifest

<strong>generateNewManifest.py</strong>

Python script which takes in a mutect results directory (from mutect pipeline found [here](https://github.com/GavinHaLab/LNCaP-subline-paper/blob/main/pipelines/mutect/mutect2_bq.snakefile)) and MANE canonical transcript file (used for sorting annovar-annotated amino-acid changes based on transcript ID) and produces a filtered SNV and Indel manifest file for use in downstream genomic analyses

<strong>getCounts.py</strong>

Script used to obtain specific variant type counts seen in tables 1, 2 from SNV/Indel manifest produced by generateNewManifest.py

<strong>getUniqueCounts.py</strong>

Separate script used to obtain the SNV/Indel overlaps seen between sublines, values used in table 2


<strong>makePhyloTable.py</strong>

Script used to generate table for phylogenetic tree seen in figure 2C. Input is SNV/Indel manifest, output is matrix with rows as variant IDs, and columns containing some variant-level information (function, annotations, etc.) as well as one col per sample, with values 1 (present) or 0 (absent) for each sample. 

<strong>getVennCountsV2.py</strong>

Script which takes in output from makePhyloTable.py, determines number of SNVs shared between samples pairs, all lines of interest (in our case, FGC lines), and unique to each line. Used to inform counts shown in S1G


