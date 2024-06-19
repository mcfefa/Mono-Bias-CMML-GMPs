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