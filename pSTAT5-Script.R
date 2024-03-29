### R Script analyzing pSTAT5 TotalSeq pilot and projecting onto BCD dataset

### DIRECTORIES

# BCD data: /blue/ferrallm/00_data/single-cell/CMML/BCD/allSeurat_39+8_postStandardPipeline_withHarmony_withSeuratGeneScores_05-20-2021_withWNT-2022-04-26.rds

# TotalSeq pSTAT5 Pilot
# Parent Directory: /blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4563-Padron
#    2-C-001_unstimulated
#    2-C-001_GM-CSF
#    KB-14_unstimulated
#    KB-14_GM-CSF

setwd("~/Mono-Bias-CMML-GMPs")

libraryPath <- "/home/ferrallm/Mono-Bias-CMML-GMPs/lib2"

## installing Seurat packages
# remotes::install_github("mojaveazure/seurat-object", "seurat5", lib=libraryPath)
# remotes::install_github("satijalab/seurat", "seurat5", lib=libraryPath, quiet = TRUE)

library('ggplot2', lib.loc=libraryPath)
library('Seurat', lib.loc=libraryPath)
library('patchwork', lib.loc=libraryPath)
library('forcats', lib.loc=libraryPath)
library('tidyverse', lib.loc=libraryPath)
library('data.table', lib.loc=libraryPath)

# remotes::install_github("satijalab/seurat-data", "seurat5", lib=libraryPath, quiet = TRUE)
# remotes::install_github("satijalab/azimuth", "seurat5", lib=libraryPath, quiet = TRUE)
# remotes::install_github("satijalab/seurat-wrappers", "seurat5", lib=libraryPath, quiet = TRUE)
# remotes::install_github("stuart-lab/signac", "seurat5", lib=libraryPath, quiet = TRUE)

## output strings and directories
dir <- "/blue/ferrallm/00_data/single-cell/CMML/totalseq-results/pSTAT5/pSTAT5-TotalSeq-Pilot-"
date <- "2024-03-28"


##################################################################
##  CREATING INDIVIDUAL RDS FILES FOR EACH TOTALSEQ PATIENT
##################################################################

savedir <- "/blue/ferrallm/00_data/single-cell/CMML/totalseq-results/pSTAT5/pSTAT5-TotalSeq-Pilot-"
enddir <-"Unprocessed-RDS_2024-02-13.rds"

##### SAMPLE 1
P1file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4563-Padron/2-C-001_unstimulated/outs/per_sample_outs/2-C-001_unstimulated/count/sample_filtered_feature_bc_matrix"

P1.counts <- Read10X(data.dir=P1file)
P1.rna <- P1.counts$`Gene Expression`
P1.adt <- P1.counts$`Antibody Capture`

all.equal(colnames(P1.rna), colnames(P1.adt))

# creates a Seurat object based on the scRNA-seq data
P1 <- CreateSeuratObject(counts=P1.rna, project="2C001unstim")

# create a new assay to store ADT information
# P1_adt_assay <- CreateAssay5Object(counts = P1.adt) <--- errored this time, unclear why
P1_adt_assay <- CreateAssayObject(counts=P1.adt)

# add this assay to the previously created Seurat object
P1[["ADT"]] <- P1_adt_assay

# Validate that the object now contains multiple assays
Assays(P1)

# Extract a list of features measured in the ADT assay
rownames(P1[["ADT"]])

P1 <- RenameCells(object=P1, add.cell.id="2C001unstim")
saveRDS(P1, paste(savedir,"2-C-001-unstimulated",enddir,sep=""))

####### SAMPLE 2
P2file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4563-Padron/2-C-001_GM-CSF/outs/per_sample_outs/2-C-001_GM-CSF/count/sample_filtered_feature_bc_matrix"
PtName <- "2C001stim"

P2.counts <- Read10X(data.dir=P2file)
P2.rna <- P2.counts$`Gene Expression`
P2.adt <- P2.counts$`Antibody Capture`

P2 <- CreateSeuratObject(counts=P2.rna, project=PtName)
P2_adt_assay <- CreateAssayObject(counts = P2.adt)
P2[["ADT"]] <- P2_adt_assay
P2 <- RenameCells(object=P2, add.cell.id=PtName)
saveRDS(P2, paste(savedir,"2-C-001-GMCSF-stimulated",enddir,sep=""))

####### SAMPLE 3
P3file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4563-Padron/KB-14_unstimulated/outs/per_sample_outs/KB-14_unstimulated/count/sample_filtered_feature_bc_matrix"
Pt3Name <- "KB14unstim"

P3.counts <- Read10X(data.dir=P3file)
P3.rna <- P3.counts$`Gene Expression`
P3.adt <- P3.counts$`Antibody Capture`

P3 <- CreateSeuratObject(counts=P3.rna, project=Pt3Name)
P3_adt_assay <- CreateAssayObject(counts = P3.adt)
P3[["ADT"]] <- P3_adt_assay
P3 <- RenameCells(object=P3, add.cell.id=Pt3Name)
saveRDS(P3, paste(savedir,"KB-14-unstimulated",enddir,sep=""))

####### SAMPLE 4
P4file <- "/blue/ferrallm/00_data/single-cell/CMML/Moffitt-CICPT-4563-Padron/KB-14_GM-CSF/outs/per_sample_outs/KB-14_GM-CSF/count/sample_filtered_feature_bc_matrix"
Pt4Name <- "KB14stim"

P4.counts <- Read10X(data.dir=P4file)
P4.rna <- P4.counts$`Gene Expression`
P4.adt <- P4.counts$`Antibody Capture`

P4 <- CreateSeuratObject(counts=P4.rna, project=Pt4Name)
P4_adt_assay <- CreateAssayObject(counts = P4.adt)
P4[["ADT"]] <- P4_adt_assay
P4 <- RenameCells(object=P4, add.cell.id=Pt4Name)
saveRDS(P4, paste(savedir,"KB14-GMCSF-stimulated",enddir,sep=""))

####### MERGE SEURAT OBJECTS
TotalSeqCohort <- merge(x=P1, y=list(P2,P3,P4), merge.data=TRUE, project="pSTAT5totalSeq")
saveRDS(TotalSeqCohort, paste(savedir,"Unprocessed-Merged",enddir,sep=""))

##################################################################
##  QUALITY CONTROL - FILTER, NORMALIZE, SCALE DATA
##################################################################

## TotalSeqCohort[['']] <- TotalSeqCohort@active.ident

######## Get Summary of nFeature RNA Distribution
summary(TotalSeqCohort@meta.data$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.0   240.0   389.0   440.1   565.0  3185.0 

######## CMML TotalSeq Stats
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
# [1] 18082 34429

## filtered based on approach
dim(subset(TotalSeqCohort, subset = nFeature_RNA > nFeatLower_CMML & nFeature_RNA < nFeatUpper_CMML & percent.mito < perMitoUpper_CMML))
# [1] 18082 12356

## filtered based on strict number values from BCD
dim(subset(TotalSeqCohort, subset = nFeature_RNA > nFeatLower_CMML & nFeature_RNA < 5808.62 & percent.mito < perMitoUpper_CMML))
# [1] 18082 13870

######## Filter Cohort based on cutoffs
### decided just to filter based on mitochondrial content because otherwise losing over half of the data and there doesn't look to be a ton of outlying cells
TotalSeqCohort <- subset(TotalSeqCohort, subset = percent.mito < perMitoUpper_CMML)
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

## Load TotalSeq Cohort Seurat Object (R crashed)
# totalseqfile <- "/blue/ferrallm/00_data/single-cell/CMML/totalseq-results/CMML-TotalSeq-Cohort-PostScaling_2023-12-15.rds"
# TotalSeqCohort <- readRDS(totalseqfile)

# Run, plot and save PCA
TotalSeqCohort <- RunPCA(TotalSeqCohort, features = VariableFeatures(object = TotalSeqCohort))

# PC_ 1 
# Positive:  MSI2, EGFL7, ITM2C, PBXIP1, PRSS57, MTURN, PTMS, NPDC1, FAM30A, PPFIBP1 
# HSPG2, SPNS2, CAVIN1, OBSCN, SNCA, ZSCAN18, DIPK1B, COL6A2, BAALC, ZFP36L2 
# CHMP6, MEIS1, FAM118A, REXO2, ARVCF, CLEC11A, SPINK2, MEGF6, BEX3, CPXM1 
# Negative:  ITGAX, LYZ, IFI30, VEGFA, TNFRSF1B, INSIG1, IER3, PPIF, IL1B, CD300E 
# OLR1, PLAUR, SOD2, SLC2A3, SRA1, APBB3, CTSB, KYNU, IL1RN, SEMA6B 
# G0S2, MMP19, PFKFB3, THBD, AQP9, VCAN, S100A9, CXCL8, OSM, LPXN 
# PC_ 2 
# Positive:  JUN, HSPA1B, RHOB, IER5, FOS, CCL5, PMAIP1, HSPA1A, IFIT2, KLF2 
# HERC5, SYNE2, GADD45G, DNAJB1, OASL, GNLY, SPOCK2, IFIH1, ID2, HSPH1 
# TNFAIP3, IRF1, BICDL1, ZAP70, IL32, SOCS1, ZC3HAV1, KLRG1, ISG15, TNIK 
# Negative:  LMNA, ANPEP, EGFL7, SEMA6B, PPIF, FAM30A, HSPG2, MTURN, BAALC, ITM2C 
# CRIP2, DPYSL3, ADA, SPNS2, PTMS, MSI2, SCARF1, CAVIN1, CHMP6, PPFIBP1 
# TNFRSF1B, PRSS57, ADGRG1, DIPK1B, AHRR, EHD2, ELN, RAB13, ARVCF, CD9 
# PC_ 3 
# Positive:  RHOB, GADD45G, IER2, JUN, HSPA1A, H1FX, IER5, ID2, IER5L, JUND 
# HSPA1B, HIST1H1B, SOX4, CITED4, FOS, MAFB, ID1, IRS2, KLF4, H2AFX 
# ASF1B, CEBPA, EGR1, RASD1, NUSAP1, SOCS1, DNAJB1, HSPA2, H1F0, AZU1 
# Negative:  BICDL1, SYNE2, HERC5, IFIT2, SPOCK2, ZC3H12D, IL7R, OASL, ZAP70, LTB 
# TNIK, CD6, IL32, ZC3HAV1, GNLY, IFIT3, ITPKB, TCF7, ITK, PLAAT4 
# ETS1, CD8A, TRAC, TNFRSF25, STAT4, TRBC1, CD2, BCL11B, KLRK1, TRBC2 
# PC_ 4 
# Positive:  HIST1H1B, HBB, HBD, HIST1H1C, MYBL2, SLC4A1, HMBS, SLC2A1, PPIF, AHSP 
# SPTB, HBM, MPO, SPTA1, THBD, HMGB2, NUSAP1, ANK1, HIST1H1E, PRC1 
# CA1, SHCBP1, RRM2, SPC24, CTSG, BTG2, CA2, OSM, MKI67, PRDX2 
# Negative:  EPPK1, HSPA1B, VPS37B, HSPH1, CCDC86, ZNF165, CRABP2, ITGAD, CDC42BPG, HSPA6 
# SLC4A1AP, PLEKHA6, HSPA1A, DCUN1D3, PPP1R13L, CTRC, ALDOC, CREM, PPFIBP1, TMEM240 
# ANKRD28, BGLAP, ERI2, NCR1, NABP1, C17orf107, MATN1, BTBD9, DNAJB6, SPAG5 
# PC_ 5 
# Positive:  HBB, HBD, SLC4A1, HMBS, AHSP, SPTB, HBM, SPTA1, ANK1, CA1 
# DMTN, RIPOR3, SELENBP1, GATA1, FCAR, PRDX2, SLC25A37, EPPK1, UBAC1, CA2 
# KEL, HBA2, BLVRB, ACHE, ALAS2, TNXB, ERI2, RXRA, ALAD, EPB42 
# Negative:  RHOB, IER2, GADD45G, JUND, JUN, ID2, IER5L, RASD1, IER5, MAFB 
# SOCS1, CEBPA, HEXIM1, NRARP, ID1, BTG2, CITED4, IRF7, EGR1, H2AFX 
# SOX4, KLF2, FOS, PHLDA2, ZFP36L2, CDKN1C, CD74, CSKMT, SOD2, GADD45B 

pdf(paste(dir, "PCA_DimPlot_", date,".pdf",sep=""), width = 11, height = 6)
DimPlot(TotalSeqCohort, reduction = "pca", pt.size = 0.0001)
dev.off()

# Determine dimensionality of dataset
pdf(paste(dir,"PCA_ElbowPlot_",date,".pdf",sep=""), width = 11, height = 6)
ElbowPlot(TotalSeqCohort, ndims = 50)
dev.off()

saveRDS(TotalSeqCohort, paste(dir,"Cohort-PostPCA_",date,".rds",sep=""))

##################################################################
## MULTIMODAL REFERENCE MAPPING
##################################################################
## Tutorial: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

## Load BCD Seurat Object
reffile <- "/blue/ferrallm/00_data/single-cell/CMML/BCD/allSeurat_39+8_postStandardPipeline_withHarmony_withSeuratGeneScores_05-20-2021_withWNT-2022-04-26.rds"
reference <- readRDS(reffile)

## Load TotalSeq Cohort Seurat Object 
# totalseqfile <- "/blue/ferrallm/00_data/single-cell/CMML/totalseq-results/pSTAT5/pSTAT5-TotalSeq-Pilot-Cohort-PostPCA_2024-02-13.rds"
# TotalSeqCohort <- readRDS(totalseqfile)

## BATCHING SAMPLES AND THEN FIND ANCHORS AND MAP
TotalSeqCohort.batches <- SplitObject(TotalSeqCohort, split.by = "orig.ident")

## Computing a cached neighbor index
reference <- FindNeighbors(
  object = reference,
  reduction = "pca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

## Find Anchors
anchors <- list()
for (i in 1:length(TotalSeqCohort.batches)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = reference,
    query = TotalSeqCohort.batches[[i]],
    k.filter = NA,
    reference.reduction = "pca",
    reference.neighbors = "spca.annoy.neighbors",
    dims = 1:50
  )
}
saveRDS(anchors, paste(dir,"Anchors-for-BCD_Individual-Batches_",date,".rds",sep=""))

# anchorfile <- "/blue/ferrallm/00_data/single-cell/CMML/totalseq-results/CMML-TotalSeq-Anchors-for-BCD_Individual-Batches_2023-12-19.rds"
# anchors <- readRDS(anchorfile)

## Individual Mapping
##   running/troubleshooting with firrst sample and then will loop through remainder of batches
##   also slimmed down to just repredicting clustering -- re-doing for a wnn.umap to align with tutorial

reference <- RunUMAP(reference, dims=1:50, reduction.name = "umap.v2", reduction.key = "UMAPv2_", return.model = TRUE)

## adding in BCD UMAP and clusters explicitly as meta-data so can hopefully visualize using these in the future
##    can also predict clusters from this meta data and confirm assignment
BCDumapfile <- "/blue/ferrallm/00_data/single-cell/CMML/BCD/BCD-UMAP-Embeddings_2020-03-31.csv"
BCD.embeddings <- read.csv(BCDumapfile)
reference[['BCD-UMAP-1']] <- BCD.embeddings$UMAP_1[match(rownames(reference@meta.data),BCD.embeddings$X)]
reference[['BCD-UMAP-2']] <- BCD.embeddings$UMAP_2[match(rownames(reference@meta.data),BCD.embeddings$X)]

BCDclustersfile <- "/blue/ferrallm/00_data/single-cell/CMML/BCD/BCD-UMAP-Clusters_2020-03-31.csv"
BCD.clusters <- read.csv(BCDclustersfile)
reference[['BCD-Clusters']] <- BCD.clusters$RNA_snn_res.0.05[match(rownames(reference@meta.data),BCD.embeddings$X)]

saveRDS(reference, paste(dir,"Reference-BCD-Updated-UMAPv2+origEmbed_",date,".rds",sep=""))

DimPlot(reference, group.by = "clusterResolution_0.05", reduction = "umap.v2") 

TotalSeqCohort.batches[[1]] <- MapQuery(
  anchorset = anchors[[1]], 
  query = TotalSeqCohort.batches[[1]],
  reference = reference, 
  refdata = list(
    predicted_cluster = "clusterResolution_0.05",
    predicted_BCD_UMAP1 = "BCD-UMAP-1",
    predicted_BCD_UMAP2 = "BCD-UMAP-2",
    predicted_BCD_Clusters = "BCD-Clusters"
    ),
  reference.reduction = "pca",
  reduction.model = "umap.v2"
)

## test visualization --- reference vs predicted
p1 <- DimPlot(reference, reduction = 'umap.v2', group.by = 'clusterResolution_0.05', label.size = 3)
p2 <- DimPlot(TotalSeqCohort.batches[[1]], reduction = 'ref.umap', group.by = 'predicted.predicted_cluster', label.size = 3)
p1 + p2 + plot_layout(guides = "collect")

DimPlot(reference, dim=vector(refUMAP[,2], refUMAP[,3]), group.by = 'clusterResolution_0.05', label.size = 3)

a <- reference$`BCD-UMAP-1`
m <- matrix(unlist(a),byrow=TRUE,ncol=length(a[[1]]))
rownames(m) <- names(a)
as.data.frame(m)

b <- reference$`BCD-UMAP-2`
n <- matrix(unlist(b),byrow=TRUE,ncol=length(b[[1]]))
rownames(n) <- names(b)
as.data.frame(n)

refUMAP <- merge(m,n,by="row.names",all=TRUE)
c <- as.data.frame(refUMAP)


for (i in 2:length(TotalSeqCohort.batches)) {
  TotalSeqCohort.batches[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = TotalSeqCohort.batches[[i]],
    reference = reference, 
    refdata = list(
      predicted_cluster = "clusterResolution_0.05"),
    reference.reduction = "pca",
    reduction.model = "umap.v2"
  )
}

saveRDS(TotalSeqCohort.batches, paste(dir,"Cohort-mapped-to-BCD-attempt2_",date,".rds",sep=""))


date <- "2024-02-15"
reference <- readRDS("/blue/ferrallm/00_data/single-cell/CMML/totalseq-results/pSTAT5/pSTAT5-TotalSeq-Pilot-Reference-BCD-Updated-UMAPv2+origEmbed_2024-02-13.rds") 
TotalSeqCohort.batches <- readRDS("/blue/ferrallm/00_data/single-cell/CMML/totalseq-results/pSTAT5/pSTAT5-TotalSeq-Pilot-Cohort-mapped-to-BCD-attempt2_2024-02-13.rds") 


## VISUALIZATION
pdf(paste(dir, "DimPlot_predicted-clustering_", date, ".pdf",sep=""), width = 18, height = 12)
p1 <- DimPlot(reference, reduction = 'umap.v2', group.by = 'clusterResolution_0.05', label.size = 3)
p2 <- DimPlot(TotalSeqCohort.batches[[1]], reduction = 'ref.umap', group.by = 'predicted.predicted_cluster', label.size = 3)
p3 <- DimPlot(TotalSeqCohort.batches[[2]], reduction = 'ref.umap', group.by = 'predicted.predicted_cluster', label.size = 3)
p4 <- DimPlot(TotalSeqCohort.batches[[3]], reduction = 'ref.umap', group.by = 'predicted.predicted_cluster', label.size = 3)
p5 <- DimPlot(TotalSeqCohort.batches[[4]], reduction = 'ref.umap', group.by = 'predicted.predicted_cluster', label.size = 3)
pTotal <- p1 + p2 + p3 + p4 + p5 + plot_layout(guides = "collect")
print(pTotal)
dev.off()

## re-merge batches
TotalSeqCohortm <- merge(TotalSeqCohort.batches[[1]], TotalSeqCohort.batches[2:length(TotalSeqCohort.batches)], merge.dr = "ref.umap")
saveRDS(TotalSeqCohortm, paste(dir,"Cohort-mapped-to-BCD-attempt2-merged_",date,".rds",sep=""))

pdf(paste(dir, "DimPlot_predicted-clustering_merged_", date, ".pdf",sep=""), width = 18, height = 12)
dPlot <- DimPlot(TotalSeqCohortm, reduction = "ref.umap", group.by =  "predicted.predicted_cluster", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
print(dPlot)
dev.off()

pdf(paste(dir, "DimPlot_BCD-clustering_UMAPv2_", date, ".pdf",sep=""), width = 18, height = 12)
dPlot <- DimPlot(reference, reduction = "umap.v2", group.by =  "clusterResolution_0.05", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
print(dPlot)
dev.off()

## exporting cluster composition to csv
file1 <- paste(dir,"outputDataNames_",date,".csv",sep="")
file2 <- paste(dir,"outputData_",date,".csv",sep="")

# save the cluster identity for each individual cell --- this is the predicted_cluster based on BCD clustering
cat(TotalSeqCohortm@meta.data$predicted.predicted_cluster, file=file2, sep=",\n")

# save the orig.ident per cell (this is the totalseq patient name)
listNames <- TotalSeqCohortm@meta.data$orig.ident
write.table(data.frame(listNames),
            row.names = FALSE,
            col.names = FALSE, 
            file = file1,
            sep = ",")

# merging the two files together in grouping to ultimately write out the number of cells per identity per cluster 
#(so in this case, we'd end up with A being a table that has 5 samples x number of cluster idenfied number of rows and two columns )
mydat1 <- read.csv(file2)
mydat2 <- read.csv(file1)
fulldat <- cbind(mydat2[1],mydat1[1])
fulltab <- as.data.table(fulldat)

# name table columns
names(fulltab)[1] <- paste("patient")
names(fulltab)[2] <- paste("cluster")

# group data based on clusters
group_by(fulltab, cluster)

# create a table counting unqiue UMIs/cells per cluster
tabPerClus <- fulltab %>% group_by(cluster) %>% count()
type <- sub("\\_.*","",fulltab$patient)
fulltab <- cbind(fulltab, type)
A <- fulltab %>% group_by(cluster) %>% count(type)

#Order A and sum same components
A <- A[order(A$type),]
count <- 0
for (i in 1:length(A$type)){
  if ( i == 1){
    A[i,'PartialSum'] <- A[i,'n']
    count <- count+1
  }else if (A[i, 'type'] %in% A[i-1, 'type']){
    A[i, 'PartialSum'] <- A[i-1, 'PartialSum'] + A[i, 'n']
    count <- count + 1
    if (i == length(A$type)){
      A[c((i-count):(i)), 'Sum'] <- A[i, 'PartialSum']
      break
    }
  }else{
    A[i, 'PartialSum'] <- A[i, 'n']
    if (i-1-count > 0){
      A[c((i-1-count):(i-1)), 'Sum'] <- A[i-1, 'PartialSum']
    }else{
      A[c((i-count):(i-1)), 'Sum'] <- A[i-1, 'PartialSum']
    }
    count <- 0
  }
}

A$Fraction <- A$n/A$Sum

# saving matrix/table A as a CSV file that will later be read into Julia for diversity calculations
divout <- paste(dir,"CellBreakdown_PerClusterPerType_predicted-cluster_",date,".csv",sep="") 
write.csv(A, file=divout)

##################################################################
## VISUALIZATION OF ADT INFORMATION
##################################################################

## Marker identified from COMET analysis with TotalSeq Panel B ADT
##    CD14, ITGAM, FCGR2A, FCGR1A, CD48, CD36
## 
##    search alternative gene names: https://www.ncbi.nlm.nih.gov/gene/

## list of antibodies
rownames(TotalSeqCohortm[["ADT"]])
# [1] "anti-STAT5"

## seems like ADT data needs some processing
TotalSeqCohortm <- NormalizeData(TotalSeqCohortm, normalization.method="CLR", margin=2, assay="ADT")
saveRDS(TotalSeqCohortm, paste(dir,"Cohort-mapped-to-BCD-attempt2-merged-postADTprocessing_",date,".rds",sep=""))

## pSTAT5
pdf(paste(dir, "FeaturePlot_pSTA5_ref-v-rna-v-adt_", date, ".pdf",sep=""), width = 18, height = 6)
p1 <- FeaturePlot(reference, features = c("STAT5A"), reduction = "umap.v2", label.size = 3)
p2 <- FeaturePlot(TotalSeqCohortm, features = c("rna_STAT5A"), reduction = "ref.umap", label.size = 3)
p3 <- FeaturePlot(TotalSeqCohortm, features = c("adt_anti-STAT5"), reduction = 'ref.umap', max.cutoff = 3)
pp <- p1 | p2 | p3 
print(pp)
dev.off()

## CD38
pdf(paste(dir, "FeaturePlot_CD38_ref-v-rna_", date, ".pdf",sep=""), width = 18, height = 6)
p1 <- FeaturePlot(reference, features = c("CD38"), reduction = "umap.v2", label.size = 3)
p2 <- FeaturePlot(TotalSeqCohortm, features = c("rna_CD38"), reduction = "ref.umap", label.size = 3)
pp <- p1 | p2
print(pp)
dev.off()

## CD34
pdf(paste(dir, "FeaturePlot_CD34_ref-v-rna_", date, ".pdf",sep=""), width = 18, height = 6)
p1 <- FeaturePlot(reference, features = c("CD34"), reduction = "umap.v2", label.size = 3)
p2 <- FeaturePlot(TotalSeqCohortm, features = c("rna_CD34"), reduction = "ref.umap", label.size = 3)
pp <- p1 | p2
print(pp)
dev.off()

## ADT 2-C-001
pdf(paste(dir, "FeaturePlot_ADT_2-C-001-unstim-v-stim_", date, ".pdf",sep=""), width = 18, height = 6)
p1 <- FeaturePlot(TotalSeqCohort.batches[[1]], features = c("adt_anti-STAT5"), reduction = "ref.umap", label.size = 3) #, title=TotalSeqCohort.batches[[2]]@active.ident)
p2 <- FeaturePlot(TotalSeqCohort.batches[[2]], features = c("adt_anti-STAT5"), reduction = "ref.umap", label.size = 3) #, title=TotalSeqCohort.batches[[1]]@active.ident)
pp <- p1 | p2
print(pp)
dev.off()

pdf(paste(dir, "VlnPlot_ADT_2-C-001-unstim-v-stim_", date, ".pdf",sep=""), width = 6, height = 18)
p1 <- VlnPlot(TotalSeqCohort.batches[[1]], features = c("adt_anti-STAT5")) 
p2 <- VlnPlot(TotalSeqCohort.batches[[2]], features = c("adt_anti-STAT5"))
pp <- p1 + p2
print(pp)
dev.off()

## ADT KB-14
pdf(paste(dir, "FeaturePlot_ADT_KB-14-unstim-v-stim_", date, ".pdf",sep=""), width = 18, height = 6)
p1 <- FeaturePlot(TotalSeqCohort.batches[[3]], features = c("adt_anti-STAT5"), reduction = "ref.umap", label.size = 3) #, title=TotalSeqCohort.batches[[2]]@active.ident)
p2 <- FeaturePlot(TotalSeqCohort.batches[[4]], features = c("adt_anti-STAT5"), reduction = "ref.umap", label.size = 3) #, title=TotalSeqCohort.batches[[1]]@active.ident)
pp <- p1 | p2
print(pp)
dev.off()

##################################################################
## pSTAT5 ANALYSIS
##################################################################

reference <- readRDS("/blue/ferrallm/00_data/single-cell/CMML/totalseq-results/pSTAT5/pSTAT5-TotalSeq-Pilot-Reference-BCD-Updated-UMAPv2+origEmbed_2024-02-13.rds") 
TotalSeqCohortm <- readRDS("/blue/ferrallm/00_data/single-cell/CMML/totalseq-results/pSTAT5/pSTAT5-TotalSeq-Pilot-Cohort-mapped-to-BCD-attempt2-merged_2024-02-15.rds") 

TotalSeqCohortm <- NormalizeData(TotalSeqCohortm, normalization.method="CLR", margin=2, assay="ADT")
saveRDS(TotalSeqCohortm, paste(dir,"Cohort-mapped-to-BCD-attempt2-merged-postADTprocessing_",date,".rds",sep=""))

## pSTAT5
pdf(paste(dir, "FeaturePlot_pSTA5_ref-v-rna-v-adt_", date, ".pdf",sep=""), width = 18, height = 6)
p1 <- FeaturePlot(reference, features = c("STAT5A"), reduction = "umap.v2", label.size = 3)
p2 <- FeaturePlot(TotalSeqCohortm, features = c("rna_STAT5A"), reduction = "ref.umap", label.size = 3)
p3 <- FeaturePlot(TotalSeqCohortm, features = c("adt_anti-STAT5"), reduction = 'ref.umap', max.cutoff = 3)
pp <- p1 | p2 | p3 
print(pp)
dev.off()

## ADT only
pdf(paste(dir, "FeaturePlot_pSTA5_adt_", date, ".pdf",sep=""), width = 6, height = 6)
p3 <- FeaturePlot(TotalSeqCohortm, features = c("adt_anti-STAT5"), reduction = 'ref.umap')
print(p3)
dev.off()

DefaultAssay(TotalSeqCohortm) <- "ADT"

### visualizing ADT based on sample -- DotPlot
pdf(paste(dir, "DotPlot_pSTA5_adt_bySample_", date, ".pdf",sep=""), width = 6, height = 6)
p4 <- DotPlot(TotalSeqCohortm, features = c("adt_anti-STAT5")) #, group.by = TotalSeqCohortm$predicted.predicted_cluster)
print(p4)
dev.off()

### visualizing ADT based on sample -- VlnPlot
pdf(paste(dir, "VlnPlot_pSTA5_adt_bySample_", date, ".pdf",sep=""), width = 6, height = 6)
p5 <- VlnPlot(TotalSeqCohortm, features = c("adt_anti-STAT5"))
print(p5)
dev.off()

### visualizing ADT based on predicted cluster --- DotPlot
pdf(paste(dir, "DotPlot_pSTA5_adt_byBCDcluster_", date, ".pdf",sep=""), width = 6, height = 6)
p4 <- DotPlot(TotalSeqCohortm, features = c("adt_anti-STAT5"), group.by = "predicted.predicted_cluster")
print(p4)
dev.off()

### visualizing ADT based on predicted cluster --- VlnPlot
pdf(paste(dir, "VlnPlot_pSTA5_adt_byBCDcluster_", date, ".pdf",sep=""), width = 6, height = 6)
p4 <- VlnPlot(TotalSeqCohortm, features = c("adt_anti-STAT5"), group.by = "predicted.predicted_cluster")
print(p4)
dev.off()

### visualizing ADT split by cluster and grouped by sample type ---- Vln Plot
pdf(paste(dir, "VlnPlot_pSTA5_adt_", date, ".pdf",sep=""), width = 10, height = 6)
p6 <- VlnPlot(TotalSeqCohortm, features = c("adt_anti-STAT5"), split.by = "predicted.predicted_cluster")
print(p6)
dev.off()

divout <- paste(dir,"Data_Sample+ClusterSplit-for-VlnPlot_",date,".csv",sep="") 
write.csv(p6[[1]]$data, file=divout)

### visualizing ADT split by cluster and grouped by sample type ---- Dot Plot
pdf(paste(dir, "DotPlot_pSTA5_adt_", date, ".pdf",sep=""), width = 10, height = 6)
p7 <- DotPlot(TotalSeqCohortm, features = c("adt_anti-STAT5"), split.by = "predicted.predicted_cluster", cols='RdBu')
print(p7)
dev.off()

divout2 <- paste(dir,"Data_Sample+ClusterSplit-for-DotPlot_",date,".csv",sep="") 
write.csv(p7$data, file=divout2)


##<--------------------------------------------------------
