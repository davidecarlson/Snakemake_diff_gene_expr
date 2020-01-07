import shutil

configfile:
    "config_hosp100319.yaml"
print(config['data'])

SAMPLES, = glob_wildcards(config['data']+"/{id}_R1.fastq.gz")
print(SAMPLES)

RESULTS = config['results_loc']
print(RESULTS)

onstart:
	print("The samples used in this analysis are:"),
	print(SAMPLES),
	print("Running the differential gene expression (DGE) pipeline") 

rule all:
	input:
		expand(config['results_loc'] + "/results/fastq_trimmed/{sample}.1.trimmed.fastq", sample=SAMPLES),
		expand(config['results_loc'] + "/results/salmon/{sample}", sample=SAMPLES),
		config['results_loc'] + "/results/deseq2/DEseq2_results.RData"

include:	"rules/preprocess.smk"
include:	"rules/salmon.smk"
include:	"rules/deseq2.smk"

onsuccess:
	print("Differential Gene Expression analysis finished!")
	shutil.rmtree(".snakemake")
