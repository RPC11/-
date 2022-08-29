install.packages("rlist")
install.packages("Seurat")
install.packages("SCENIC")
library(rlist)
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
library(tidyverse)
library(data.table)
library(RcisTarget)
library(SCENIC)

set.seed(221)
##
#
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
#######################################################
######################################################
### Co-expression network
#test <- pmbc.updated@assays$RNA@counts
#scenicOptions
#org <- "hgnc" #mgi  hgnc
#dbDir <- "cisTarget_databases" # RcisTarget databases location
#myDatasetTitle <- "SCENIC example on Tcell" # choose a name for your analysis
#data(defaultDbNames)
#dbs <- defaultDbNames[[org]]
########################################
library(RcisTarget)
#download the feather files
#dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather",
 #            "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs
#dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
#for(featherURL in dbFiles)
#{
 # download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
#}
library(SCENIC)
scenicOptions <- initializeScenic(org= "hgnc" , dbDir="cisTarget_databases" , nCores=5) 
exprMat<-immune.combined@assays$RNA@counts


############################################################################## 

#Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)






# Optional: Binarize activity
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Export:
# saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)


#################################################################
############################################################################
########################################################################
#Expression matrix of the top five genes from the top five cells
immune.combined@assays$RNA@counts[1:5,1:5]
#Cells and genes after PCA analysis were visualized
VizDimLoadings(immune.combined, dims = 1:2, reduction = "pca")
#head(pbmc@reductions$pca@cell.embeddings)
#immune.combined[["RNA"]]@counts

























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







