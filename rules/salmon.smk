
if config['species'] == "human":
	INDEX="/gpfs/projects/GenomicsCore/indexes/human_gencode_v36/salmon_index"
if config['species'] == "mouse":
	INDEX="/gpfs/projects/GenomicsCore/indexes/mouse_gencode_m25/salmon_index"
if config['species'] == "neither":
	INDEX=config['index']

#INDEX = config['index']

if config['paired_reads'] == "True":

	rule salmon_pe:
		input:
			fwd=config['results_loc']+"/results/fastq_trimmed/{sample}.1.trimmed.fastq",
			rev=config['results_loc']+"/results/fastq_trimmed/{sample}.2.trimmed.fastq",
			index=INDEX
		output:
			quant=directory(config['results_loc']+"/results/salmon/{sample}")
		threads: int(config['threads']) 
		priority:   50
		shell: "salmon quant -i {input.index} -l A -1 {input.fwd} -2 {input.rev} --validateMappings --threads {threads} --seqBias --gcBias -o {output.quant}"

if config['paired_reads'] == "False":
	rule salmon_se:
		input:
			reads=config['results_loc']+"/results/fastq_trimmed/{sample}.trimmed.fastq",
			index=INDEX
		output:
			quant=directory(config['results_loc']+"/results/salmon/{sample}")
		threads: int(config['threads'])
		priority:   50
		shell:	"salmon quant -i {input.index} -l A -r {input.reads} --validateMappings --threads {threads} --seqBias --gcBias -o {output.quant}"
