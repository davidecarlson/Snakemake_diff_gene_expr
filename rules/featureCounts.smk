
if config['species'] == "human":
	ANNOTATION="/gpfs/projects/GenomicsCore/indexes/human_gencode_v36/gencode.v36.annotation.gtf"
if config['species'] == "mouse":
	ANNOTATION="/gpfs/projects/GenomicsCore/indexes/mouse_gencode_m25/gencode.vM25.annotation.gtf"
if config['species']=="neither":
	ANNOTATION=config['gtf']

#FASTQs=config['my_fastqs']

rule featureCounts:
	input:
		bams=expand(config['results_loc'] + "/results/STAR/passTwo/{sample}/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
	output:
		countsfile=config['results_loc'] + "/results/featureCounts/counts.txt",
		condensed=config['results_loc'] + "/results/featureCounts/condensed_counts.txt",
		lengths=config['results_loc'] + "/results/featureCounts/gene_lengths.txt"
	threads: 20
	params:
		gtf=ANNOTATION
	log:
		config['results_loc'] + "/results/logs/featureCounts.log"
		
	shell: "featureCounts -a {params.gtf} -t exon -g gene_id -o {output.countsfile} -T {threads} {input.bams} 2> {log} "
		"&& cut -f1,6- {output.countsfile} > {output.condensed} "
		"&& cat {output.condensed} | sort -n -k1,1 | cut -f2 | head -n -2 > {output.lengths}"
