#manta Normal Tumor Paired snakefile
#Michael Yang and Akira Nair
#Template created 06/2023
#Ha Lab
#Fred Hutchinson Cancer Research Center

# Snakefile to run manta mutation caller
# https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md
# ^ User Guide ^
# See Run Configuration & Execution section for details
"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Java/11.0.2
ml Python

#command to run snakemake (remove -np at end when done validating):
snakemake -s manta_tumorOnly.snakefile --latency-wait 120 --restart-times 0 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 60 -np

"""
configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input:
        expand("results_122023/{tumors}/runWorkflow.py", tumors=config["samples"]),
        expand("results_122023/{tumors}/results/variants/tumorSV.vcf.gz", tumors=config["samples"]),
        expand("results_122023/{tumors}/results/variants/candidateSV.vcf.gz", tumors=config["samples"]),
        expand("results_122023/{tumors}/results/variants/candidateSmallIndels.vcf.gz", tumors=config["samples"]),
        expand("results_122023/{tumors}/results/stats/alignmentStatsSummary.txt", tumors=config["samples"]),
        expand("results_122023/{tumors}/results/stats/svLocusGraphStats.tsv", tumors=config["samples"]),
        expand("results_122023/{tumors}/results/stats/svCandidateGenerationStats.tsv", tumors=config["samples"]),
        expand("results_122023/{tumors}/results/stats/svCandidateGenerationStats.xml", tumors=config["samples"])

def split_sample(wildcards):
	return expand("{sample}", sample=(((config["samples"][wildcards.tumors][0]).split('/')[-1])).split('.')[0])


# manta configuration
rule manta_config:
    input:
        tumor_filepath = lambda wildcards: config["samples"][wildcards.tumors]
        #normal_filepath = lambda wildcards: config["samples"][wildcards.tumors][2]
    output:
        manta_config = protected("results_122023/{tumors}/runWorkflow.py")
    params:
        reference_fasta = config["reference_fasta"],
        manta = config["manta_install"],
        interpreter = config["interpreter"],
        #tumor_name = lambda wildcards: config["samples"][wildcards.tumors],
        project_dir = expand(config["project_dir"]),
        tumor_name = "results_122023/{tumors}"
    log:
        "logs/manta_config/{tumors}.manta_config.log"
    shell:
        "({params.interpreter} {params.manta}/bin/configManta.py \
        --tumorBam {input.tumor_filepath} \
        --referenceFasta {params.reference_fasta} \
        --runDir {params.project_dir}/{params.tumor_name}/) 2> {log}"


# manta execution
rule manta_run:
    input:
        manta_workflow="results_122023/{tumors}/runWorkflow.py"
    output:
        ("results_122023/{tumors}/results/variants/tumorSV.vcf.gz"),
        ("results_122023/{tumors}/results/variants/candidateSV.vcf.gz"),
        ("results_122023/{tumors}/results/variants/candidateSmallIndels.vcf.gz"),
        ("results_122023/{tumors}/results/stats/alignmentStatsSummary.txt"),
        ("results_122023/{tumors}/results/stats/svLocusGraphStats.tsv"),
        ("results_122023/{tumors}/results/stats/svCandidateGenerationStats.tsv"),
        ("results_122023/{tumors}/results/stats/svCandidateGenerationStats.xml")
    params:
        interpreter = config["interpreter"]
    log:
        "logs/manta_run/{tumors}.manta_run.log"
    shell:
        "({params.interpreter} {input.manta_workflow} -j 10 ) 2> {log}"
