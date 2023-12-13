### R Script for projecting TotalSeq patient samples onto BCD dataset

### DIRECTORIES

# BCD data: /blue/ferrallm/00_data/single-cell/CMML/BCD/allSeurat_39+8_postStandardPipeline_withHarmony_withSeuratGeneScores_05-20-2021_withWNT-2022-04-26.rds

# TotalSeq 
# Batch1: /blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-01
#    P1_RX_1_001_CMML_MPN_NRAS
#    P2_RX_4_001_CMML_MPN_KRAS
#    P3_RX_8_001_CMML_MPN_KRAS
#    P4_KB_13_103_010_CMML_NRAS
#    P5_KB_14_103_011_CMML_NRAS
# Batch2: /blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-02
#    P6_1_X_001_CMML_MPN_KRAS
#    P7_2_U_001_CMML_MPN_NRAS
#    P8_KB_1_003_101_CMML_MPN_NRAS
#    P9_KB_10_103_007_CMML_NRAS
#    P10_PDX_9_002_CMML_MDS_KRAS

setwd("~/Mono-Bias-CMML-GMPs")

libraryPath <- "/home/ferrallm/Mono-Bias-CMML-GMPs/lib2"

## installing Seurat packages
remotes::install_github("mojaveazure/seurat-object", "seurat5", lib=libraryPath)
remotes::install_github("satijalab/seurat", "seurat5", lib=libraryPath, quiet = TRUE)

library('Seurat', lib.loc=libraryPath)


remotes::install_github("satijalab/seurat-data", "seurat5", lib=libraryPath, quiet = TRUE)
remotes::install_github("satijalab/azimuth", "seurat5", lib=libraryPath, quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", lib=libraryPath, quiet = TRUE)
remotes::install_github("stuart-lab/signac", "seurat5", lib=libraryPath, quiet = TRUE)


#### CREATING INDIVIDUAL RDS FILES FOR EACH TOTALSEQ PATIENT

savedir <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-01/CMML-"
enddir <-"_Unprocessed-RDS_2023-12-13.rds"

P1file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-01/P1_RX_1_001_CMML_MPN_NRAS/outs/filtered_feature_bc_matrix/"

P1.counts <- Read10X(data.dir=P1file)
P1.rna <- P1.counts$`Gene Expression`
P1.adt <- P1.counts$`Antibody Capture`

all.equal(colnames(P1.rna), colnames(P1.adt))

##<------------------------------
# creates a Seurat object based on the scRNA-seq data
P1 <- CreateSeuratObject(counts=P1.rna, project="RX1001")

# create a new assay to store ADT information
P1_adt_assay <- CreateAssay5Object(counts = P1.adt)

# add this assay to the previously created Seurat object
P1[["ADT"]] <- P1_adt_assay

# Validate that the object now contains multiple assays
Assays(P1)

# Extract a list of features measured in the ADT assay
rownames(P1[["ADT"]])


##class(P1[["RNA"]])
P1 <- RenameCells(object=P1, add.cell.id="RX1001")
saveRDS(P1, paste(savedir,"P1_RX-1-001_CMML-MPN-NRAS",enddir,sep=""))

