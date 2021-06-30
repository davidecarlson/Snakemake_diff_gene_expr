rule fpkm:
	input:
		script="/gpfs/projects/GenomicsCore/Snakemake_diff_gene_expr/scripts/fpkm.R",
		counts=expand(config['results_loc'] + "/results/htseq/{sample}_rawcounts.csv", sample=SAMPLES),
		lengths=config['results_loc']+"/results/featureCounts/gene_lengths.txt",
	output:
		fpkm=config['results_loc']+"/results/deseq2/fpkm.csv",

	log:	config['results_loc']+"/results/logs/deseq2.log"
	params:
		counts_dir=config['results_loc']+"/results/htseq/",
		outdir=config['results_loc']+"/results/deseq2/"
	priority:	0
	shell:
		"Rscript {input.script} {params.counts_dir} {input.lengths} {params.outdir} 2> {log}"
