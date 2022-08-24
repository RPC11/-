install.packages("rlist")
install.packages("Seurat")
library(rlist)
library(Matrix)

library(Seurat)
library(tidyverse)
library(data.table)
library(Seurat)
library(tidyverse)
library(data.table)

set.seed(221)
##
#
raw_count_list <- readRDS("E:/work/Tcell/activated_T_cells/activated_T_cells.rds")
metadata <- read.csv("E:/work/Tcell/activated_T_cells/metadata_activated_T_cells.csv")

count_combined <- rlist::list.cbind(raw_count_list)
#Initialize the Seurat object
state_seurat <- CreateSeuratObject(count_combined, project = "SeuratProject", assay = "RNA",
                                   min.cells = 10, min.features = 1, meta.data = NULL)
# Initialize the Seurat object
state_seurat[["percent.mt"]] <- PercentageFeatureSet(state_seurat, pattern = "^MT-")

#data has been QCed
##The number of genes, the number of UMI, and the proportion of mitochondrial genes were plotted

pdf(file= "QC_metrics.pdf", width = 20, height = 11)
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


# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
#Run nonlinear dimensionality reduction
#Seurat provides several nonlinear dimensionality reduction techniques, 
#such as TSNE and UMAP, to visualize and explore these datasets. 
#We further reduced the data in the PC dimension of the previous step to two dimensions to visualize the data. 
#Previously, PCA was linear dimensionality reduction, which could deal with high-dimensional data well. 
#Tsne and UMAP were nonlinear dimensionality reduction methods, which were fast for low-dimensional data, 
#but difficult to deal with high-dimensional data. Compared with TSNE, 
#UMAP can better reflect the original data structure and has better continuity.
pdf(file = "E:/work/Tcell/UMAP_activated_T_State.pdf", width = 25, height = 18)
h<-DimPlot(immune.combined, group.by = "State")
h
dev.off()
pdf(file = "E:/work/Tcell/UMAP_activated_T_Tissue.pdf", width = 25, height = 18)
h1<-DimPlot(immune.combined, group.by = "Tissue")
h1
dev.off()



class(immune.combined@assays$RNA@counts)

############################################################################

#Expression matrix of the top five genes from the top five cells
immune.combined@assays$RNA@counts[1:5,1:5]
#Cells and genes after PCA analysis were visualized
VizDimLoadings(immune.combined, dims = 1:2, reduction = "pca")
#head(pbmc@reductions$pca@cell.embeddings)
#immune.combined[["RNA"]]@counts

###############################################
install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes =TRUE)
####################################
head(immune.combined[["RNA"]]@dcounts)
BiocManager::install("SCopeLoomR")
install.packages("SCopeLoomR")
library(SCopeLoomR)
build_loom(file.name = "tcell.loom",
           dgem = immune.combined@assays$RNA@counts,
           title = "tcellloom",
           default.embedding = immune.combined@reductions$umap.rna@cell.embeddings,
           default.embedding.name = "umap.rna")



























############################################
#The installation of SCENIC
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
# If your bioconductor version is previous to 4.0, see the section bellow

## Required
BiocManager::install(c("AUCell", "RcisTarget"),ask = F,update = F) 
BiocManager::install(c("GENIE3"),ask = F,update = F)  # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"),ask = F,update = F) 
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"),ask = F,update = F)
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"),ask = F,update = F)
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")
###############################################


library(GENIE3)

head(immune.combined[["RNA"]]@counts)
immune.combined[["RNA"]]@counts[1:5,1:5]
dim(immune.combined[["RNA"]]@counts)

#########################################
network = grnboost2(expression_data=ex_matrix,
                    tf_names=tf_names)



#########################
#GENIE3  Transcriptional regulatory network analysis
set.seed(123)
weightMat <- GENIE3(as.matrix(immune.combined[["RNA"]]@counts[1:50,1:50]),nCores=4,nTrees=50)
dim(weightMat)
weightMat
#############################################
library(SCopeLoomR)
##Parameters are sufficient to add information：
loom <- open_loom(file.name)
add_hierarchy(loom = loom, hierarchy = create_hierarchy(level.1.name = "Mouse", level.2.name = "Toy Datasets", level.3.name = ""))
add_col_attr(loom=loom, key = "Cell type", value=cell.info$cellType, as.annotation=T)
###Add Seurat information
seurat.annotation<-read.table(file = paste0(seuratDir, "Res2_Clusters.tsv", header=T, quote = '', sep = "\t", stringsAsFactors=F))
add_seurat_clustering(loom = loom
                      , seurat = seurat
                      , default.clustering.resolution = "res.2"
                      , annotation = seurat.annotation
                      , annotation.cluster.id.cn = "res.2" 
                      , annotation.cluster.description.cn = "Annotation")


