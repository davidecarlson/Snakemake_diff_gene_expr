library(stringr)
library(readr)
library(assertr)
library(tximport)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Mm.eg.db)
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
gtf <- args[4]
samples_input <- args[5]

#Make a transcript database object use the gencode35 GTF file. GTF file location specified based on config file
txdb <- makeTxDbFromGFF(gtf, "gtf", "gencode35")

#Prepare a tx2gene dataframe to associate transcript names with gene names
k <- keys(txdb, keytype="TXNAME")	
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Read in sample-specific information from a previously generated table.  Currently must create this manually
samples <- read.table(samples_input, header = T)
#get list of salmon quant files and confirm they exist
dir <- getwd()
files <-file.path(indir, samples$Run, "quant.sf")
names(files) <- paste0("sample_", samples$Run)
all(file.exists(files))

#get abundances for each sample
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar=TRUE)

# prepare for DESeq2
ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples,design = ~ Condition) # This design should be modified as needed to fit different experimental setups

# get fpkm normalized results for each gene

fpkm_res <- as.data.frame(fpkm(ddsTxi, robust = TRUE))

# convert Ensembl IDs to gene name and symbols for fpkm results

ens.str_fpkm <- substr(rownames(fpkm_res), 1, 18)

print(ens.str_fpkm)
print(typeof(fpkm_res))

fpkm_res$geneName <- mapIds(org.Mm.eg.db,
                     keys=ens.str_fpkm,
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

fpkm_res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ens.str_fpkm,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

#print(fpkm_res)

# save fpkm results to output file

write.csv(as.data.frame(fpkm_res),file=paste(outdir, "fpkm.csv",sep=""))

# filter rows were total counts are less than 10
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

# add gene name and symbol to the results

ens.str <- substr(rownames(res), 1, 18)

#print(ens.str)

res$geneName <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")


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

#report_pdf <- DESeq2Report(ddsTxi, project = 'DESeq2 PDF report', 
#    intgroup = c('Condition'), outdir = outdir,
#    output = 'DESeq2_Report', theme = theme_bw(), output_format = 'pdf_document',
#    device = 'pdf')
    
# save results to csv file
write.csv(as.data.frame(resOrdered),file=paste(outdir, "DEseq2_results.csv",sep=""))

## Save the workspace
save.image(file=paste(outdir, "DEseq2_results.RData",sep=""))
