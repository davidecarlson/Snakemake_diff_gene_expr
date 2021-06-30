rule make_fastq_list:
	input:
		readdir=config['results_loc']+"/results/fastq_trimmed",
		#reads=config['results_loc']+"/results/fastq_trimmed/{sample}.trimmed.fastq",
	output:
		outfile=config['results_loc']+"/results/fastq_trimmed/fastq_list.tsv"
	shell:	"python Snakemake_diff_gene_expr/scripts/make_fastq_list.py --inPath {input.readdir} --output {output.outfile}"
