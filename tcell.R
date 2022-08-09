install.packages("rlist")
library(rlist)
library(Matrix)
install.packages("Seurat")
library(Seurat)
library(tidyverse)
library(data.table)


set.seed(221)

#import activated_T_cells.rds and metadata_activated_T_cells.csv
files_list<-readRDS("E:/work/Tcell/activated_T_cells/activated_T_cells.rds")
metadata<-read.csv("E:/work/Tcell/activated_T_cells/metadata_activated_T_cells.csv")
count_combined<-rlist::list.cbind(files_list)

#Create count_combined data as a Seurat object
state_seurat <- CreateSeuratObject(count_combined, project = "SeuratProject", assay = "RNA",
                                   min.cells = 10, min.features = 1, meta.data = NULL)
state_seurat[["percent.mt"]] <- PercentageFeatureSet(state_seurat, pattern = "^MT-")

#data has been QCed, we get the state
pdf(file= "QC_metrics.pdf", width = 20, height = 11)
p <- VlnPlot(state_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p
dev.off()

# Add metadata
seurat_meta <- state_seurat@meta.data
#Assign the line name of Meta. Data in state_seurate to Cellid
seurat_meta$cellid <- rownames(seurat_meta)

metadata$orig.ident <- gsub("_.*","",metadata$Patient)
seurat_meta <- left_join(seurat_meta, metadata, by = "orig.ident")
rownames(seurat_meta) <- seurat_meta$cellid
state_seurat@meta.data <- seurat_meta

# To do from here

# split the dataset into a list of two seurat objects (stim and CTRL)
seurat_list <- SplitObject(state_seurat, split.by = "Tissue")

# normalize and identify variable features for each dataset independently
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
#Rank the features
features <- SelectIntegrationFeatures(object.list = seurat_list)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#After homogenization of single sample, multiple samples were integrated
immune.anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, reduction = "rpca")
immune.combined <- IntegrateData(anchorset = immune.anchors)


# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
export
pdf(file = "E:/work/Tcell/activated_T_cells/UMAP_activated_T_State.pdf", width = 25, height = 18)
DimPlot(immune.combined, group.by = "State")
dev.off()
pdf(file = "E:/work/Tcell/activated_T_cells/UMAP_activated_T_Tissue.pdf", width = 25, height = 18)
DimPlot(immune.combined, group.by = "Tissue")
dev.off()