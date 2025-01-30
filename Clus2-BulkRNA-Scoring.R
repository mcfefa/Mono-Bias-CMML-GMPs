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
#install.packages("BiocManager", lib=libraryPath)
library(BiocManager)#, lib.loc=libraryPath)
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

# BiocManger::install("GSVA")
library(GSVA)
## bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html

setwd("/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/CITESeq_BulkRNASeq_Results_05312024")

## read in genelists
glfile <- "/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2+HSC-GeneSignatures.csv"
genelists <- as.list(read.csv(glfile,sep=","))

## run GSVA
clus2_es <- gsva(bulkcts, genelists)

library(RColorBrewer)

pdf(file="/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2_GSVA-enrichment-score-hm_2025-01-30.pdf", width=10, height=10)
heatmap(clus2_es, 
        margins=c(15,15))
dev.off()

## restrict to Clus2, HSC, GMP, MEP
genelists2 <- genelists
genelists2[5] <- NULL # removed Wu ProB
genelists2[5] <- NULL # removed Wu ETP

## run GSVA
clus2_es2 <- gsva(bulkcts, genelists2)

pdf(file="/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2_GSVA-enrichment-score-HSC-GMP-MEP-hm_2025-01-30.pdf", width=10, height=10)
heatmap(clus2_es2, 
        margins=c(15,15))
dev.off()

## loading in genesets from Surendra
library(GSEABase)
genesetAtt2 <- getGmt("/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/combined_Clus2_BCD_human_hematopoiesis.gmt")

## run GSVA
clus2_es2 <- gsva(bulkcts, genesetAtt2)

library(RColorBrewer)

pdf(file="/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2_GSVA-enrichment-score-hm-bloodGenelists_2025-01-30.pdf", width=10, height=10)
heatmap(clus2_es2, 
        margins=c(15,15))
dev.off()

write.csv(clus2_es2, "/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2_GSVA-enrichment-score-hm-scores-bloodGenelists_2025-01-30.csv")

### updated with different version of refined Clus2 signatures
##    subset from the original "clus2.de.markers.csv" for upregulated (>0 FC cutoffs and adj p-val <0.05) -- this was from Clus2 versus all other cells
##    clus2all1426: all genes with FC > 0.25, which gives a signature with 1426 genes
##    clus2pt5FC741: all genes with FC > 0.5, which gives a signature with 741 genes
##    clus2w1FC276: all genes with FC > 1, which gives a signature with 276 genes
##    clus2w1p5FC126: all genes with FC > 1.5, which gives a signature with 126 genes
##    clus2w2FC63: all genes with FC > 2, which gives a signature with 63 genes
##    clus2w2p5FC34: all genes with FC > 2.5, which gives a signature with 34 genes

glfileO <- "/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2-GeneSignatures-Opt_2025-01-30.csv"
genelistsO <- as.list(read.csv(glfileO,sep=","))

clus2_esO <- gsva(bulkcts, genelistsO)
pdf(file="/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2_GSVA-enrichment-score-HSC-GMP-MEP-hm_OptClus2_2025-01-30.pdf", width=10, height=10)
heatmap(clus2_esO, 
        margins=c(15,15))
dev.off()

write.csv(clus2_esO, "/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2_GSVA-enrichment-score-hm-scores-Optclus2_2025-01-30.csv")

### update with a verison of refined Clus2 signatures based on DGE of Clus2 vs Clus0
##    subset from original "clus2-vs-clus0.de.markers.csv" in the revision directory for upregulated (>0 FC cutoffs and adj p-val <0.05)
##    clus2v0all1490: all genes with FC > 0.25, which gives a signature with 1490 genes
##    clus2v0pt5FC767: all genes with FC > 0.5, which gives a signature with 767 genes
##    clus2v0w1FC288: all genes with FC > 1, which gives a signature with 288 genes
##    clus2v0w1p5FC126: all genes with FC > 1.5, which gives a signature with 126 genes
##    clus2v0w2FC64: all genes with FC > 2, which gives a signature with 64 genes
##    clus2v0w2p5FC37: all genes with FC > 2.5, which gives a signature with 37 genes

glfileOatt2 <- "/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2-GeneSignatures-Opt2v0_2025-01-30.csv"
genelistsOatt2 <- as.list(read.csv(glfileOatt2,sep=","))

clus2_esOatt2 <- gsva(bulkcts, genelistsOatt2)
pdf(file="/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2_GSVA-enrichment-score-HSC-GMP-MEP-hm_OptClus2v0_2025-01-30.pdf", width=10, height=10)
heatmap(clus2_esOatt2, 
        margins=c(15,15))
dev.off()

write.csv(clus2_esOatt2, "/blue/ferrallm/00_data/RNAseq/Moffitt-CICPT_4448_Padron_RNAseq-TRE/Clus2_GSVA-enrichment-score-hm-scores-Optclus2v0_2025-01-30.csv")

### future steps --- 
## refine GMT files to include only the significant 'other' signatures and our BCD signatures for visualization
gs1 <- GeneSetCollection(genelistsO$clus2w2FC63)
updGeneSet <- GeneSetCollection()















