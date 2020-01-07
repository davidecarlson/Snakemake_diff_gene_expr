INDEX = config['index']

rule salmon:
	input:
		fwd=config['results_loc']+"/results/fastq_trimmed/{sample}.1.trimmed.fastq",
		rev=config['results_loc']+"/results/fastq_trimmed/{sample}.2.trimmed.fastq",
		index=INDEX
	output:
		quant=directory(config['results_loc']+"/results/salmon/{sample}")
	threads: 28
	shell: "salmon quant -i {input.index} -l A -1 {input.fwd} -2 {input.rev} --validateMappings --threads {threads} --seqBias --gcBias -o {output.quant}"
