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

##################################################################
##  CREATING INDIVIDUAL RDS FILES FOR EACH TOTALSEQ PATIENT
##################################################################

savedir <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-01/CMML-"
enddir <-"_Unprocessed-RDS_2023-12-15.rds"

##### PATIENT 1
P1file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-01/P1_RX_1_001_CMML_MPN_NRAS/outs/filtered_feature_bc_matrix/"

P1.counts <- Read10X(data.dir=P1file)
P1.rna <- P1.counts$`Gene Expression`
P1.adt <- P1.counts$`Antibody Capture`

all.equal(colnames(P1.rna), colnames(P1.adt))

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

P1 <- RenameCells(object=P1, add.cell.id="RX1001")
saveRDS(P1, paste(savedir,"P1_RX-1-001_CMML-MPN-NRAS",enddir,sep=""))

####### PATIENT 2
P2file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-01/P2_RX_4_001_CMML_MPN_KRAS/outs/filtered_feature_bc_matrix/"
PtName <- "RX4001"

P2.counts <- Read10X(data.dir=P2file)
P2.rna <- P2.counts$`Gene Expression`
P2.adt <- P2.counts$`Antibody Capture`

P2 <- CreateSeuratObject(counts=P2.rna, project=PtName)
P2_adt_assay <- CreateAssay5Object(counts = P2.adt)
P2[["ADT"]] <- P2_adt_assay
P2 <- RenameCells(object=P2, add.cell.id=PtName)
saveRDS(P2, paste(savedir,"P2_RX-4-001_CMML-MPN-KRAS",enddir,sep=""))

####### PATIENT 3
P3file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-01/P3_RX_8_001_CMML_MPN_KRAS/outs/filtered_feature_bc_matrix/"
Pt3Name <- "RX8001"

P3.counts <- Read10X(data.dir=P3file)
P3.rna <- P3.counts$`Gene Expression`
P3.adt <- P3.counts$`Antibody Capture`

P3 <- CreateSeuratObject(counts=P3.rna, project=Pt3Name)
P3_adt_assay <- CreateAssay5Object(counts = P3.adt)
P3[["ADT"]] <- P3_adt_assay
P3 <- RenameCells(object=P3, add.cell.id=Pt3Name)
saveRDS(P3, paste(savedir,"P3_RX-8-001_CMML-MPN-KRAS",enddir,sep=""))

####### PATIENT 4
P4file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-01/P4_KB_13_103_010_CMML_NRAS/outs/filtered_feature_bc_matrix/"
Pt4Name <- "KB13103010"

P4.counts <- Read10X(data.dir=P4file)
P4.rna <- P4.counts$`Gene Expression`
P4.adt <- P4.counts$`Antibody Capture`

P4 <- CreateSeuratObject(counts=P4.rna, project=Pt4Name)
P4_adt_assay <- CreateAssay5Object(counts = P4.adt)
P4[["ADT"]] <- P4_adt_assay
P4 <- RenameCells(object=P4, add.cell.id=Pt4Name)
saveRDS(P4, paste(savedir,"P4_KB-13-103-010_CMML-MPN-NRAS",enddir,sep=""))

####### PATIENT 5
P5file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-01/P5_KB_14_103_011_CMML_NRAS/outs/filtered_feature_bc_matrix/"
Pt5Name <- "KB14103011"

P5.counts <- Read10X(data.dir=P5file)
P5.rna <- P5.counts$`Gene Expression`
P5.adt <- P5.counts$`Antibody Capture`

P5 <- CreateSeuratObject(counts=P5.rna, project=Pt5Name)
P5_adt_assay <- CreateAssay5Object(counts = P5.adt)
P5[["ADT"]] <- P5_adt_assay
P5 <- RenameCells(object=P5, add.cell.id=Pt5Name)
saveRDS(P5, paste(savedir,"P5_KB-14-103-011_CMML-MPN-NRAS",enddir,sep=""))

####### PATIENT 6
P6file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-02/P6_1_X_001_CMML_MPN_KRAS/outs/filtered_feature_bc_matrix/"
Pt6Name <- "1X001"

P6.counts <- Read10X(data.dir=P6file)
P6.rna <- P6.counts$`Gene Expression`
P6.adt <- P6.counts$`Antibody Capture`

P6 <- CreateSeuratObject(counts=P6.rna, project=Pt6Name)
P6_adt_assay <- CreateAssay5Object(counts = P6.adt)
P6[["ADT"]] <- P6_adt_assay
P6 <- RenameCells(object=P6, add.cell.id=Pt6Name)
saveRDS(P6, paste(savedir,"P6_1-X-001_CMML-MPN-KRAS",enddir,sep=""))

####### PATIENT 7
P7file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-02/P7_2_U_001_CMML_MPN_NRAS/outs/filtered_feature_bc_matrix/"
Pt7Name <- "2U001"

P7.counts <- Read10X(data.dir=P7file)
P7.rna <- P7.counts$`Gene Expression`
P7.adt <- P7.counts$`Antibody Capture`

P7 <- CreateSeuratObject(counts=P7.rna, project=Pt7Name)
P7_adt_assay <- CreateAssay5Object(counts = P7.adt)
P7[["ADT"]] <- P7_adt_assay
P7 <- RenameCells(object=P7, add.cell.id=Pt7Name)
saveRDS(P7, paste(savedir,"P7_2-U-001_CMML-MPN-NRAS",enddir,sep=""))

####### PATIENT 8
P8file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-02/P8_KB_1_003_101_CMML_MPN_NRAS/outs/filtered_feature_bc_matrix/"
Pt8Name <- "KB1003101"

P8.counts <- Read10X(data.dir=P8file)
P8.rna <- P8.counts$`Gene Expression`
P8.adt <- P8.counts$`Antibody Capture`

P8 <- CreateSeuratObject(counts=P8.rna, project=Pt8Name)
P8_adt_assay <- CreateAssay5Object(counts = P8.adt)
P8[["ADT"]] <- P8_adt_assay
P8 <- RenameCells(object=P8, add.cell.id=Pt8Name)
saveRDS(P8, paste(savedir,"P8_KB-1-003-101_CMML-MPN-NRAS",enddir,sep=""))

####### PATIENT 9
P9file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-02/P9_KB_10_103_007_CMML_NRAS/outs/filtered_feature_bc_matrix/"
Pt9Name <- "KB10103007"

P9.counts <- Read10X(data.dir=P9file)
P9.rna <- P9.counts$`Gene Expression`
P9.adt <- P9.counts$`Antibody Capture`

P9 <- CreateSeuratObject(counts=P9.rna, project=Pt9Name)
P9_adt_assay <- CreateAssay5Object(counts = P9.adt)
P9[["ADT"]] <- P9_adt_assay
P9 <- RenameCells(object=P9, add.cell.id=Pt9Name)
saveRDS(P9, paste(savedir,"P9_KB-10-103-007_CMML-NRAS",enddir,sep=""))

####### PATIENT 10
P10file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4448-TotalSeq-Batch-02/P10_PDX_9_002_CMML_MDS_KRAS/outs/filtered_feature_bc_matrix/"
Pt10Name <- "PDX9002"

P10.counts <- Read10X(data.dir=P10file)
P10.rna <- P10.counts$`Gene Expression`
P10.adt <- P10.counts$`Antibody Capture`

P10 <- CreateSeuratObject(counts=P10.rna, project=Pt10Name)
P10_adt_assay <- CreateAssay5Object(counts = P10.adt)
P10[["ADT"]] <- P10_adt_assay
P10 <- RenameCells(object=P10, add.cell.id=Pt10Name)
saveRDS(P10, paste(savedir,"P10_PDX-9-002_CMML-MDS-KRAS",enddir,sep=""))

####### MERGE SEURAT OBJECTS
TotalSeqCohort <- merge(x=P1, y=list(P2,P3,P4,P5,P6,P7,P8,P9,P10), merge.data=TRUE, project="TotalSeq")
saveRDS(TotalSeqCohort, paste(savedir,"TotalSeqCohort-Merged",enddir,sep=""))

##################################################################
##  QUALITY CONTROL - FILTER, NORMALIZE, SCALE DATA
##################################################################

## TotalSeqCohort[['']] <- TotalSeqCohort@active.ident

######## Get Summary of nFeature RNA Distribution
summary(TotalSeqCohort@meta.data$nFeature_RNA)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   47    2511    3707    3698    4871   10798 

######## Original Stats from CMML Cohort published in BCD
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17    1067    2418    2496    3654    9742

## Determining how upper cutoff for nFeatures compares
## nFeatUpper_CMML <- mean(first39@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(first39@meta.data$nFeature_RNA, na.rm=TRUE) # 5808.62 (consistent with Meghan)
nFeatUpper_CMML <- mean(TotalSeqCohort@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(TotalSeqCohort@meta.data$nFeature_RNA, na.rm=TRUE) # [1] 7239.94
## moved ahead with same procedure, may come back and limit upper level of nFeatures in the future
nFeatLower_CMML <- 450 

## Calculate mito percentage
TotalSeqCohort[["percent.mito"]] <- PercentageFeatureSet(TotalSeqCohort, pattern = "^MT-")
perMitoUpper_CMML <- 25

## Visualize QC metrics
dir <- "/blue/ferrallm/00_data/single-cell/CMML/totalseq-results/CMML-TotalSeq-"
date <- "2023-12-15"

pdf(paste(dir,"VlnPlot_nFeature+nCount+PercentMito_",date,".pdf",sep=""), width = 11, height = 6)
v <- VlnPlot(TotalSeqCohort, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
print(v)
dev.off()

pdf(paste(dir, "FeatureScatter_nCountvMito_nCountvnFeature_", date, ".pdf",sep=""), width = 11, height = 6)
plot1 <- FeatureScatter(TotalSeqCohort, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(TotalSeqCohort, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cPlot <- CombinePlots(plots = list(plot1, plot2))
print(cPlot)
dev.off()

## original dataset dimensions
dim(TotalSeqCohort)
# [1]  36601 101828

## filtered based on approach
dim(subset(TotalSeqCohort, subset = nFeature_RNA > nFeatLower_CMML & nFeature_RNA < nFeatUpper_CMML & percent.mito < perMitoUpper_CMML))
# [1] 36601 92415

## filtered based on strict number values from BCD
dim(subset(TotalSeqCohort, subset = nFeature_RNA > nFeatLower_CMML & nFeature_RNA < 5808.62 & percent.mito < perMitoUpper_CMML))
# [1] 36601 82922

######## Filter Cohort based on cutoffs
TotalSeqCohort <- subset(TotalSeqCohort, subset = nFeature_RNA > nFeatLower_CMML & nFeature_RNA < nFeatUpper_CMML & percent.mito < perMitoUpper_CMML)
saveRDS(TotalSeqCohort, paste(dir,"Cohort-PostFiltering_",date,".rds",sep=""))

######## Normalizing data
TotalSeqCohort <- NormalizeData(TotalSeqCohort, normalization.method = "LogNormalize", scale.factor = 10000)

######## Identification of highly variable features
TotalSeqCohort <- FindVariableFeatures(TotalSeqCohort, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(TotalSeqCohort), 10)

# plot variable features with and without labels
pdf(paste(dir, "VariableFeaturePlot_", date, ".pdf",sep=""), width = 11, height = 6)
plot1 <- VariableFeaturePlot(TotalSeqCohort)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
cPlot2 <- CombinePlots(plots = list(plot1, plot2))
print(cPlot2)
dev.off()

saveRDS(TotalSeqCohort, paste(dir,"Cohort-PostNormalization_",date,".rds",sep=""))

######## Scale data
# Scale data (only use HVG), regressing out effects of nCountRNA and percent.mito
TotalSeqCohort <- ScaleData(TotalSeqCohort, features = VariableFeatures(TotalSeqCohort), vars.to.regress = c("nCount_RNA","percent.mito"))

saveRDS(TotalSeqCohort, paste(dir,"Cohort-PostScaling_",date,".rds",sep=""))

##################################################################
## DIMENSION REDUCTION - PCA
##################################################################

# Run, plot and save PCA
TotalSeqCohort <- RunPCA(TotalSeqCohort, features = VariableFeatures(object = TotalSeqCohort))

pdf(paste(dir, 'PCA_standard_seurat_pipeline_preHarmony_', date, '.pdf'), width = 10, height = 2)
DimPlot(TotalSeqCohort, reduction = "pca", pt.size = 0.0001)
dev.off()

# Determine dimensionality of dataset
pdf(paste(dir, 'PCA_ElbowPlot_standard_seurat_pipeline_preHarmony_', date, '.pdf'), width = 10, height = 2)
ElbowPlot(TotalSeqCohort, ndims = 50)
dev.off()

saveRDS(TotalSeqCohort, paste(dir,"Cohort-PostPCA_",date,".rds",sep=""))

##################################################################
## MULTIMODAL REFERENCE MAPPING
##################################################################
## Tutorial: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
##<--------------------------------------------------------
## Load BCD Seurat Object
reffile <- "/blue/ferrallm/00_data/single-cell/CMML/BCD/allSeurat_39+8_postStandardPipeline_withHarmony_withSeuratGeneScores_05-20-2021_withWNT-2022-04-26.rds"
reference <- readRDS(reffile)

## Find Anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = TotalSeqCohort,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

####### TODO: UPDATE REFDATA WITH PREDICTIONS TO MAKE ---- MOSTLY CLUSTER, PHENOTYPE, ETC

## Map onto Reference
TotalSeqCohort <- MapQuery(
  anchorset = anchors,
  query = TotalSeqCohort,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

##################################################################
## 
##################################################################






