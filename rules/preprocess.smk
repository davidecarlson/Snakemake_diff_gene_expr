if config['gzipped_reads'] == "True":
	gz=".gz"
if config['gzipped_reads'] == "False":
	gz=""

rule fastp:
	input:
		fwd=config['data']+"/{sample}_R1.fastq"+gz,
		rev=config['data']+"/{sample}_R2.fastq"+gz
	output:
		fwd=config['results_loc']+"/results/fastq_trimmed/{sample}.1.trimmed.fastq",
		rev=config['results_loc']+"/results/fastq_trimmed/{sample}.2.trimmed.fastq",
		html=config['results_loc']+"/results/fastq_trimmed/{sample}.html",
		json=config['results_loc']+"/results/fastq_trimmed/{sample}.json"
	threads: int(config['threads'])
	log:	config['results_loc']+"/results/logs/{sample}/preprocess.log"
	shell:
		"fastp --in1 {input.fwd} --in2 {input.rev} "
		"--out1 {output.fwd} --out2 {output.rev} --thread {threads} --cut_tail --html {output.html} --json {output.json} 2> {log}"
