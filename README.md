# Differential Gene Expression with Snakemake


*This Snakemake workflow uses Fastp to QC the reads, Salmon for transcript quantification, Tximport to convert from transcript abundance to gene abundance, and DEseq2 for differential gene expression.* 

### INSTALLATION
* Clone the repository:
  * `git clone https://github.com/davidecarlson/Snakemake_diff_gene_expr.git`
* Use conda to install create a new environment with all the necessary software from the provided yaml file:
  * `conda env create -f environment.yml`
  
### SETUP
  
  * Edit the example_config.yaml file to provide information regarding your data and experiment:
    * path to the RNA-Seq reads
    * path to save the results
    * information regarding the experimental conditions and samples being used
    * number of threads to use
    * species species ("human", "mouse", or "neither")
    * path to your Salmon index, which must be built ahead of time
        * see [here](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode) for instructions on building the index
    * path to gtf file of gene annotations, if "neither" was chosen as your species
    * path to tab-separated file with list of samples and associated experimental condtions (e.g., "samples.txt")
    * name of the condition from "samples.txt" that should act as the control
    * whether or not input reads are gzipped
    * whether or not input reads are paired
    
  * Edit the samples.txt file to include a row for each sample, including sample name, run name (can be the same as sample name if there is only one run per sample), and experimental condition for that sample
  
   
### USAGE
  * `conda activate snakemake_DGE`
  * `snakemake --cores <# of threads * # of processes to run currently> -s <path to snakefile>`

  *All output files and logs will be saved in the directory specificed in your config file*
    
### TO-DO
  * Add more automation so that fewer files need to be manually edited prior to usage
