
if config['species'] == "human":
	INDEX="/gpfs/projects/GenomicsCore/indexes/human_gencode_v36/STAR_index"
if config['species'] == "mouse":
	INDEX="/gpfs/projects/GenomicsCore/indexes/mouse_gencode_m25/STAR_index"
if config['species'] == "neither":
	INDEX=config['index']


#FASTQs=config['my_fastqs']

rule star_firstPass:
	input:
		reads=config['results_loc']+"/results/fastq_trimmed/{sample}.trimmed.fastq",
		index=INDEX
	output:
		outconfirm=config['results_loc']+"/results/STAR/{sample}.Aligned.sortedByCoord.out.bam",
		sj=config['results_loc'] + "/results/STAR/{sample}.SJ.out.tab"
	threads: int(config['threads'])
	params:
		bamtype="BAM SortedByCoordinate",
		outprefix=config['results_loc']+"/results/STAR/{sample}."
	log:
		config['results_loc']+"/results/logs/{sample}/STAR_pass1.log"
		
	shell: "STAR --readFilesIn {input.reads} --genomeDir {input.index} --runThreadN {threads} --outFileNamePrefix {params.outprefix} --outSAMtype {params.bamtype} 2> {log}"
