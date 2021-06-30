library(stringr)
library(readr)
library(assertr)
#library(tximport)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(DESeq2)
library(ggplot2)
library(regionReport)


# specify path to htseq output files

args <- commandArgs(trailingOnly = TRUE)
print(args)
directory <- args[1]
feature_lengths <- args[2]
outdir <- args[3]

#directory <- '/gpfs/projects/GenomicsCore/results/Pawan/AKO_v_WT/results/htseq'


# gene lenghts file from featureCounts has been presorted to match htseq order of gene names

lengths <- read.table(feature_lengths, header = FALSE)

#unlist(basepairs, use.names = FALSE)

basepairs <- as.numeric(as.matrix(lengths))

# create dataframe for all htseq files
sampleFiles <- list.files(directory, pattern ="*rawcounts.csv")
print(sampleFiles)
splitfiles <- strsplit(sampleFiles,"_")
sampleNames <- unlist(splitfiles)[2*(1:length(sampleFiles))-1]
print(sampleNames)

sampleCondition <- "All"
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition) 

print(sampleTable)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ 1)
#mcols(ddsHTSeq) <- DataFrame(mcols(ddsHTSeq), basepairs)

mcols(ddsHTSeq)$basepairs <- basepairs

print(ddsHTSeq)

# get fpkm normalized results for each gene

fpkm_res <- as.data.frame(fpkm(ddsHTSeq, robust = TRUE))


# convert Ensembl IDs to gene name and symbols for fpkm results

ens.str_fpkm <- substr(rownames(fpkm_res), 1, 18)

#print(ens.str_fpkm)
#print(typeof(fpkm_res))

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

print(fpkm_res)

# save fpkm results to output file

write.csv(as.data.frame(fpkm_res),file=paste(outdir, "fpkm.csv",sep=""))

