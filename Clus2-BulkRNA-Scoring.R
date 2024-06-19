## Clus2 Scoring for Bulk RNAseq analysis

## drawing inspiration from these turorials: 
## bioinformatics-core-shared-training.github.io/Bulk_RNAseq_Course_June23
## bioinformatics-core-shared-training.github.io/Bulk_RNAseq_Course_June23/Bulk_RNAseq_Course_Base/Markdowns/05_Data_Exploration.html

setwd("~/Mono-Bias-CMML-GMPs")

libraryPath <- "/home/ferrallm/Mono-Bias-CMML-GMPs/lib2"

## libraries --- attempting with default versions on HiPerGator
library(DESeq2)
library(tidyverse)
library(tximport)
library(BiocManger)
# BiocManager::install("pasilla",lib.loc=libraryPath)
# library('pasilla', lib.loc=libraryPath)

## code
datadir <- "/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/CITESeq_BulkRNASeq_Results_05312024/Salmon_Amy_CITESeq_BulkRNAseq/Filtered/"


# Import data with pasilla following along with: 
# bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DECSeq2.html#count-matrix-input

file <- paste(datadir,"Salmon_TPM_Filtered_GeneID_Data_without_Zeros.csv",sep="")
#bulkdata <- tximport(file, type="salmon", countsFromAbundance="scaledTPM", geneIDCol="GeneID")
## still troubleshooting reading in the file; waiting for more info from Surendra on scripts used to run salmon

bulkcts <- as.matrix(read.csv(file,sep=",",row.names="GeneID"))
head(bulkcts,2)

## drop TPM_ in front of each sample name
colnames(bulkcts) <- sub("TPM_","",colnames(bulkcts))
head(bulkcts,4)

BiocManger::install("GSVA")
library(GSVA)
## bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html

setwd("/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/CITESeq_BulkRNASeq_Results_05312024")

## read in genelists
glfile <- "/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2+HSC-GeneSignatures.csv"
genelists <- as.list(read.csv(glfile,sep=","))

## run GSVA
clus2_es <- gsva(bulkcts, genelists)

library(RColorBrewer)

pdf(file="/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2_GSVA-enrichment-score-hm_2024-06-19.pdf", width=10, height=10)
heatmap(clus2_es, 
        margins=c(15,15))
dev.off()

## restrict to Clus2, HSC, GMP, MEP
genelists2 <- genelists
genelists2[5] <- NULL # removed Wu ProB
genelists2[5] <- NULL # removed Wu ETP

## run GSVA
clus2_es2 <- gsva(bulkcts, genelists2)

pdf(file="/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2_GSVA-enrichment-score-HSC-GMP-MEP-hm_2024-06-19.pdf", width=10, height=10)
heatmap(clus2_es2, 
        margins=c(15,15))
dev.off()
