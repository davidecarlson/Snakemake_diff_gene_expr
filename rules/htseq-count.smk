
if config['species'] == "human":
	ANNOTATION="/gpfs/projects/GenomicsCore/indexes/human_gencode_v36/gencode.v36.annotation.gtf"
if config['species'] == "mouse":
	ANNOTATION="/gpfs/projects/GenomicsCore/indexes/mouse_gencode_m25/gencode.vM25.annotation.gtf"
if config['species']=="neither":
	ANNOTATION=config['gtf']

#FASTQs=config['my_fastqs']

rule htseq_count:
	input:
		bam=config['results_loc']+"/results/STAR/passTwo/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
	output:
		countsfile=config['results_loc']+"/results/htseq/{sample}_rawcounts.csv",
	threads: int(config['threads'])
	params:
		gtf=ANNOTATION
	log:
		config['results_loc']+"/results/logs/{sample}/htseq-count.log"
	threads: 1
		
	shell: "htseq-count -m union -s no -t exon -i gene_id -r pos -f bam {input.bam} {params.gtf} &> {log}"
		" && grep ^ENS {log} > {output.countsfile} "
