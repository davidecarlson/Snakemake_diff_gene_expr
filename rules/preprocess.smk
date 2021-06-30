if config['gzipped_reads'] == "True":
    gz=".gz"
if config['gzipped_reads'] == "False":
    gz=""

if config['paired_reads'] == "True":
        SAMPLES, = glob_wildcards(config['data']+"/{id}_R1.fastq"+gz)
if config['paired_reads'] == "False":
        SAMPLES, = glob_wildcards(config['data']+"/{id}.fastq"+gz)

if config['paired_reads'] == "True":
	rule fastp_pe:
		input:
        		fwd=config['data']+"/{sample}_R1.fastq"+gz,
        		rev=config['data']+"/{sample}_R2.fastq"+gz
		output:
        		fwd=config['results_loc']+"/results/fastq_trimmed/{sample}.1.trimmed.fastq",
        		rev=config['results_loc']+"/results/fastq_trimmed/{sample}.2.trimmed.fastq",
        		html=config['results_loc']+"/results/fastq_trimmed/{sample}.html",
        		json=config['results_loc']+"/results/fastq_trimmed/{sample}.json"
		threads:	int(config['threads'])
		log:	config['results_loc']+"/results/logs/{sample}/preprocess.log"
		priority:	50
		shell:
        		"fastp --in1 {input.fwd} --in2 {input.rev} "
        		"--out1 {output.fwd} --out2 {output.rev} --thread {threads} --cut_tail --html {output.html} --json {output.json} 2> {log}"

if config['paired_reads'] == "False":
	rule fastp_se:
		input:
			reads=config['data'] + "/{sample}.fastq"+gz
		output:
			reads=config['results_loc'] + "/results/fastq_trimmed/{sample}.trimmed.fastq",
			html=config['results_loc'] + "/results/fastq_trimmed/{sample}.html",
			json=config['results_loc'] + "/results/fastq_trimmed/{sample}.json"
		threads:	int(config['threads'])
		log:    config['results_loc'] + "/results/logs/{sample}/preprocess.log"
		shell:
			"fastp --in1 {input.reads} --out1 {output.reads} "
			"--thread {threads} --cut_tail --html {output.html} --json {output.json} 2> {log}"
