library(stringr)
library(readr)
library(assertr)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(regionReport)

# get arguments for input and output directories
args <- commandArgs(trailingOnly = TRUE)
print(args)
indir <- args[1]
outdir <- args[2]
#conditions <- args[3]
control <- args[3]

#Make a transcript database object use the gencode31 GTF file. You will need to update this for your GTF file
txdb <- makeTxDbFromGFF("/datahome/snakemake_workflows/data/GTFs/mus_musculus/gencode.vM22.annotation.gtf", "gtf", "gencode31")

#Prepare a tx2gene dataframe to associate transcript names with gene names
k <- keys(txdb, keytype="TXNAME")	
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Read in sample-specific information from a previously generated table.  Currently must create this manually
samples <- read.table("samples.txt", header = T)
#get list of salmon quant files and confirm they exist
dir <- getwd()
files <-file.path(indir, samples$Run, "quant.sf")
names(files) <- paste0("sample_", samples$Run)
all(file.exists(files))

#get abundances for each sample
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar=TRUE)

# prepare for DESeq2
ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Condition) # This design should be modified as needed to fit different experimental setups
keep <- rowSums(counts(ddsTxi)) >= 10
ddsTxi <- ddsTxi[keep,]
# specify which condition is the reference level or control
ddsTxi$Condition <- relevel(ddsTxi$Condition, ref = control)

# run DESeq2 to find differentially expressed genes
ddsTxi <- DESeq(ddsTxi)
dds_names <- resultsNames(ddsTxi)
print(dds_names)
#get results - P-values automatically FDR adjusted with 0.1 alpha
res <- results(ddsTxi)
#reorder the results by smallest to largest p-value
resOrdered <- res[order(res$pvalue),]

#Shrink the log fold change values using apeglm
resLFC <- lfcShrink(ddsTxi, coef=dds_names[2], type="apeglm")

#generate plots
pdf(file=paste(outdir,"deseq2_results.pdf",sep=""))

# Plot counts for gene with lowest fold change (down regulated gene)
plotCounts(ddsTxi, gene=which.min(res$log2FoldChange), intgroup="Condition")

# Plot counts for gene with lowest adjusted p-value (statistically significant gene)
plotCounts(ddsTxi, gene=which.min(res$padj), intgroup="Condition")

# PCA for samples
plotPCA(rlog(ddsTxi), intgroup="Condition")+theme_bw()

# Plot effect sizes using the shrunken LFC
plotMA(resLFC, ylim=c(-2,2))

# Close the graphics device
dev.off()

# Create final report
        
report <- DESeq2Report(ddsTxi, project = 'DESeq2 HTML report',
    intgroup = c('Condition'), outdir = outdir,
    output = 'DESeq2_Report', theme = theme_bw())
    
# save results to csv file
write.csv(as.data.frame(resOrdered),file=paste(outdir, "DEseq2_results.csv",sep=""))

## Save the workspace
save.image(file=paste(outdir, "DEseq2_results.RData",sep=""))
