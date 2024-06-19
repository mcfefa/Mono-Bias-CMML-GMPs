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

## code
datadir <- "/blue/ferrallm/00_data/Moffitt-CICPT_4448_Padron_RNAseq-TRE/CITESeq_BulkRNASeq_Results_05312024/Salmon_Amy_CITESeq_BulkRNAseq/Filtered/"

file <- paste(datadir,"Salmon_TPM_Filtered_GeneID_Data_without_Zeros.csv",sep="")
bulkdata <- tximport(file, type="salmon", countsFromAbundance="scaledTPM", geneIDCol="GeneID")
## still troubleshooting reading in the file; waiting for more info from Surendra on scripts used to run salmon