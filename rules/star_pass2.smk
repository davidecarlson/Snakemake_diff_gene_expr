
if config['species'] == "human":
	INDEX="/gpfs/projects/GenomicsCore/indexes/human_gencode_v36/STAR_index"
if config['species'] == "mouse":
	INDEX="/gpfs/projects/GenomicsCore/indexes/mouse_gencode_m25/STAR_index"
if config['species'] == "neither":
	INDEX=config['index']



rule star_secondPass:
	input:
		reads=config['results_loc'] + "/results/fastq_trimmed/{sample}.trimmed.fastq",
		index=INDEX,
		sj=config['results_loc'] + "/results/STAR/{sample}.SJ.out.tab"
	output:
		outconfirm=config['results_loc'] + "/results/STAR/passTwo/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	threads: int(config['threads'])
	params:
		bamtype="BAM SortedByCoordinate",
		outprefix=config['results_loc'] + "/results/STAR/passTwo/{sample}/{sample}.",
		quant="GeneCounts"
	log:
		config['results_loc'] + "/results/logs/{sample}/STAR_pass2.log"
		
	shell: "STAR --readFilesIn {input.reads} --genomeDir {input.index} --sjdbFileChrStartEnd {input.sj} --runThreadN {threads} "
		"--outFileNamePrefix {params.outprefix} --outSAMtype {params.bamtype} --quantMode {params.quant} 2> {log}"
