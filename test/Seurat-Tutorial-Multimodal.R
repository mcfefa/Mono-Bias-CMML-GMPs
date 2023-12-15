#### Demo/Tutorial of MultiModal Data with Seurat
####    Following along with : https://satijalab.org/seurat/articles/seurat5_multimodal_vignette

setwd("~/Mono-Bias-CMML-GMPs/test")

libraryPath <- "/home/ferrallm/Mono-Bias-CMML-GMPs/lib2"

library('Seurat', lib.loc=libraryPath)
options(Seurat.object.assay.version = "v5")


##################################################################
## Load in the data
##################################################################

# Load in the RNA UMI matrix

# Note that this dataset also contains ~5% of mouse cells, which we can use as negative
# controls for the protein measurements. For this reason, the gene expression matrix has
# HUMAN_ or MOUSE_ appended to the beginning of each gene.
cbmc.rna <- as.sparse(read.csv(file = "./data/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", 
                               sep = ",", header = TRUE, row.names = 1))

# To make life a bit easier going forward, we're going to discard all but the top 100 most
# highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Load in the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file = "./data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",",
                               header = TRUE, row.names = 1))

# Note that since measurements were made in the same cells, the two matrices have identical
# column names
all.equal(colnames(cbmc.rna), colnames(cbmc.adt))

##################################################################
## Setup a Seurat object, add the RNA and protein data
##################################################################

# creates a Seurat object based on the scRNA-seq data
cbmc <- CreateSeuratObject(counts = cbmc.rna)

# We can see that by default, the cbmc object contains an assay storing RNA measurement
Assays(cbmc)

# create a new assay to store ADT information
adt_assay <- CreateAssay5Object(counts = cbmc.adt)

# add this assay to the previously created Seurat object
cbmc[["ADT"]] <- adt_assay

# Validate that the object now contains multiple assays
Assays(cbmc)

# Extract a list of features measured in the ADT assay
rownames(cbmc[["ADT"]])

# Note that we can easily switch back and forth between the two assays to specify the default
# for visualization and analysis

# List the current default assay
DefaultAssay(cbmc)

# Switch the default to ADT
DefaultAssay(cbmc) <- "ADT"

DefaultAssay(cbmc)

##################################################################
## Mapping human peripheral blood cells
##################################################################
# Tutorial 2: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html#example-2-mapping-human-bone-marrow-cells
options(SeuratData.repo.use = "http://seurat.nygenome.org")

library(SeuratData,lib.loc=libraryPath)

# load reference data
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
# load query data
InstallData('hcabm40k')
hcabm40k <- LoadData(ds = "hcabm40k")

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
              reduction.key = "wnnUMAP_", return.model = TRUE)
DimPlot(bm, group.by = "celltype.l2", reduction = "wnn.umap") 

# computing an sPCA transformation
bm <- ScaleData(bm, assay = 'RNA')
bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')

# computing a cached neighbor index
bm <- FindNeighbors(
  object = bm,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)


# query dataset preprocessing
hcabm40k.batches <- SplitObject(hcabm40k, split.by = "orig.ident")

hcabm40k.batches <- lapply(X = hcabm40k.batches, FUN = NormalizeData, verbose = FALSE)

# mapping
anchors <- list()
for (i in 1:length(hcabm40k.batches)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = bm,
    query = hcabm40k.batches[[i]],
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}

# mapping individually
for (i in 1:length(hcabm40k.batches)) {
  hcabm40k.batches[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = hcabm40k.batches[[i]],
    reference = bm, 
    refdata = list(
      celltype = "celltype.l2", 
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
}
##<--------------------------------------------------------
# explore the mapping results
p1 <- DimPlot(hcabm40k.batches[[1]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3)
p2 <- DimPlot(hcabm40k.batches[[2]], reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 3)
p1 + p2 + plot_layout(guides = "collect")

# Merge the batches 
hcabm40k <- merge(hcabm40k.batches[[1]], hcabm40k.batches[2:length(hcabm40k.batches)], merge.dr = "ref.umap")
DimPlot(hcabm40k, reduction = "ref.umap", group.by =  "predicted.celltype", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

p3 <- FeaturePlot(hcabm40k, features = c("rna_TRDC", "rna_MPO", "rna_AVP"), reduction = 'ref.umap', 
                  max.cutoff = 3, ncol = 3)

# cell type prediction scores
DefaultAssay(hcabm40k) <- 'prediction.score.celltype'
p4 <- FeaturePlot(hcabm40k, features = c("CD16 Mono", "HSC", "Prog-RBC"), ncol = 3, 
                  cols = c("lightgrey", "darkred"))

# imputed protein levels
DefaultAssay(hcabm40k) <- 'predicted_ADT'
p5 <- FeaturePlot(hcabm40k, features = c("CD45RA", "CD16", "CD161"), reduction = 'ref.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 3)
p3 / p4 / p5


##################################################################
## 
##################################################################





##################################################################
## 
##################################################################




##################################################################
## 
##################################################################


