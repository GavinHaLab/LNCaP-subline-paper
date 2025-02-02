Configurations and scripts associated with analysis pipelines for this manuscript.

## Data Processing Pipelines:

### alignmentAndMetrics:
    Because some sequencing data was received as fastq files and some as aligned bams, two processing pipelines were used for the LNCaP subline project
    For the fastq files, the following pipeline was used:
    https://github.com/GavinHaLab/fastq_to_bam_paired_snakemake

    For the aligned bams, the following pipeline was used:
    https://github.com/GavinHaLab/realign_bam_paired_snakemake

    specific paramters used can be found her in alignmentAndMetrics/config.yaml

### mutect
    Mutect2 was run in tumor-only mode to call somatic mutations (SNVs and Indels) in the LNCaP sublines, and was drawn from the following pipeline:
    https://github.com/GavinHaLab/mutect2_snakemake/tree/master

    However, modifications were made to the pipeline to suit the needs of the LNCaP subline project. The modifications were:
    - relaxed filtering parameters for mutect2, particularly with regard to f1r2-min-bq, base-quality-score-threshold, and min-base-quality-score
    - switching to tumor-only mode

    Changes to the pipeline can be found in mutect/mutect2_bq.snakefile

### TitanCNA_and_SvABA_tumorOnly:
    TitanCNA and SvABA were called as part of one processing pipeline. 

    This pipeline was used to call CNAs in the LNCaP sublines, and was drawn from the following pipeline:
    https://github.com/GavinHaLab/TitanCNA_SV_WGS_tumorOnly/tree/master

    Modifications to better handle the cell-line data were made to the pipeline, which can be seen here:
    https://github.com/GavinHaLab/TitanCNA_SV_WGS/commit/487b9b46c0674cbcd0b0fca6892ffc480157bdeb

    Specific parameters that were changed can be found in TitanCNA_and_SvABA_tumorOnly/config_hg38.yaml

    Additionally, the svaba.snakefile at https://github.com/GavinHaLab/TitanCNA_SV_WGS_tumorOnly/blob/master/svaba.snakefile was used to call structural variants in the LNCaP sublines. 
    
    Specific parameters used can be found in config_hg38.yaml
    

### structuralVariants:
<strong> Gridss2 </strong>:
    Gridss2 was run in tumor-only mode using the following pipeline:
    https://github.com/GavinHaLab/gridss2_somatic_sv_call_tumor_only

    Specific parameters can be found in structuralVariants/gridss2/config.yaml

<strong> Manta </strong>:
    Manta was run in tumor only mode using structuralVariants/manta/manta_tumorOnly.snakefile

    Specific parameters and inputs can be found in structuralVariants/manta/config.yaml

