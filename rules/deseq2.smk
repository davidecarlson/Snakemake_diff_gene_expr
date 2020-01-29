#modify conditions and control for your specific experiment

rule deseq2:
    input:
        script="scripts/runDESeq2.R",
        files=expand(config['results_loc']+"/results/salmon/{sample}",sample=SAMPLES)
        #check=config['results_loc']+"/results/salmon/{sample}.ok"
        #check=config['results_loc']+"/results/salmon/{sample}.salmon.ok"
    output:
        config['results_loc']+"/results/deseq2/DEseq2_results.RData",
    log:    config['results_loc']+"/results/logs/deseq2.log"
    params:
        indir=config['results_loc']+"/results/salmon",
        outdir=config['results_loc']+"/results/deseq2/",
        conditions="Condition_Duct3_vs_Tuft2",
        control="Tuft2"
    priority:   0
    shell:
        "echo {input.files} && Rscript {input.script} {params.indir} {params.outdir} {params.conditions} {params.control} 2> {log}"
