install.packages("rlist")
install.packages("Seurat")
install.packages("SCENIC")
library(rlist)
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
library(tidyverse)
library(RcisTarget)
library(SCENIC)

set.seed(221)
##
raw_count_list <- readRDS("E:/work/Tcell/activated_T_cells/activated_T_cells.rds")
metadata <- read.csv("E:/work/Tcell/activated_T_cells/metadata_activated_T_cells.csv")

count_combined <- rlist::list.cbind(raw_count_list)
#Initialize the Seurat object
state_seurat <- CreateSeuratObject(count_combined, project = "SeuratProject", assay = "RNA",
                                   min.cells = 10, min.features = 1, meta.data = NULL)
#############
state_seurat@assays

# Initialize the Seurat object
state_seurat[["percent.mt"]] <- PercentageFeatureSet(state_seurat, pattern = "^MT-")
head(state_seurat[["percent.mt"]])
#data has been QCed
##The number of genes, the number of UMI, and the proportion of mitochondrial genes were plotted
#pdf(file= "QC_metrics.pdf", width = 20, height = 11)
p <- VlnPlot(state_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p
dev.off()
# Add metadata
seurat_meta <- state_seurat@meta.data
seurat_meta$cellid <- rownames(seurat_meta)

metadata$orig.ident <- gsub("_.*","",metadata$Patient)
seurat_meta <- left_join(seurat_meta, metadata, by = "orig.ident")
rownames(seurat_meta) <- seurat_meta$cellid
state_seurat@meta.data <- seurat_meta

# To do from here
seurat_blood <- subset(state_seurat, subset = Tissue %in% "Blood")
seurat_list <- SplitObject(seurat_blood, split.by = "Patient")
# split the dataset into a list of two seurat objects (stim and CTRL)
#seurat_list <- SplitObject(state_seurat, split.by = "Tissue")

# normalize and identify variable features for each dataset independently
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = seurat_list)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, reduction = "rpca")
immune.combined <- IntegrateData(anchorset = immune.anchors)

#Setting it to integrated can be interpreted as: This is the consolidated data. You can think of it as batch adjusted data.
#2. Set DefaultAssay to "RNA", which means that the following analysis will be based on the original value
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
#。
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
immune.combined<-RunTSNE(immune.combined, reduction = "pca", dims = 1:30)
###########################################################
h<-DimPlot(immune.combined, group.by = "State")
h
#Use TSNE for classification
t<-DimPlot(immune.combined,reduction = "tsne",group.by = "State")
t
# Plot a legend to map colors to expression levels 
FeaturePlot(immune.combined, features = c("CD8A","CD8B", "CD4"))



##################################################################################
####################################################################################3

# The 10 genes with the most dramatic changes in expression
top10 <- head(VariableFeatures(immune.combined), 10) 
top10
#Look at the gene sets that have a high impact on each principal component
print(immune.combined[["pca"]], dims = 1:5, nfeatures = 5)
#pca
DimPlot(immune.combined, reduction = "pca",split.by = 'ident')
#heatmap
DimHeatmap(immune.combined, dims = 1, cells = 500, balanced = TRUE)


# find markers for every cluster compared to all remaining cells, report only the positive ones
immune.combined.markers  <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune.combined.markers   %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

#pbmc.markers    immune.combined.markers    
top10 <-  immune.combined.markers   %>% group_by(cluster) %>% top_n(n = 14, wt = avg_log2FC)
DoHeatmap(immune.combined, features = top10$gene) + NoLegend()

new.cluster.ids <- c("IL7R",	"IFIT3",	"LTA",	"SRM",	"BNIP3",	"HSPE1",
                     "PRF1",	"YBX3",	"PTMS",	"IFI6",	"EGR1",	"CCL19",	"HLA−DPA1", "GZMA",	"ACTB")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#DEG
