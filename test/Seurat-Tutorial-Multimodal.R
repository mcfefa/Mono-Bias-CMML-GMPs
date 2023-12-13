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
## 
##################################################################
##<--------------------------------------------------------




##################################################################
## 
##################################################################





##################################################################
## 
##################################################################




##################################################################
## 
##################################################################


