rule add_names:
	input:
		script="/gpfs/projects/GenomicsCore/Snakemake_diff_gene_expr/scripts/geneid2name.py",
		results=config['results_loc']+"/results/deseq2/DEseq2_results.csv"
	output:
		results=config['results_loc']+"/results/deseq2/DEseq2_results_with_names.tsv"
	log:
		config['results_loc']+"/results/logs/addnames.log"
	threads:	40
	shell:
		"python {input.script} --input {input.results} --threads {threads} --output {output.results} 2> {log}"
