#modify control for your specific experiment

if config['species'] == "human":
	ANNOTATION="/gpfs/projects/GenomicsCore/indexes/human_gencode_v36/gencode.v36.annotation.gtf"
if config['species'] == "mouse":
	ANNOTATION="/gpfs/projects/GenomicsCore/indexes/mouse_gencode_m25/gencode.vM25.annotation.gtf"
if config['species']=="neither":
	ANNOTATION=config['gtf']

MYSAMPLES=config['my_samples']

rule deseq2:
	input:
		script="/gpfs/projects/GenomicsCore/Snakemake_diff_gene_expr/scripts/runDESeq2.R",
		files=expand(config['results_loc']+"/results/salmon/{sample}",sample=SAMPLES)
	output:
		config['results_loc']+"/results/deseq2/DEseq2_results.RData",
	log:	config['results_loc']+"/results/logs/deseq2.log"
	params:
		indir=config['results_loc']+"/results/salmon",
		outdir=config['results_loc']+"/results/deseq2/",
		gtf=ANNOTATION,
		in_samples=MYSAMPLES,
       		control=config['control']
	priority:	0
	shell:
		"echo {input.files} && Rscript {input.script} {params.indir} {params.outdir} {params.control} {params.gtf} {params.in_samples} 2> {log}"
