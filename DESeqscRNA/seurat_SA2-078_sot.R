# *************************************************
# Seurat script. started by Subia Ahmed on February 28, 2024 
# **************************************************

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Load the datasets
SA2_078S.data <- Read10X(data.dir = "/projects/b1042/GoyalLab/Subia/20240202_scRNAseq/run_count_SA2_078_SotA/outs/filtered_feature_bc_matrix/")

#plot directory
plotDirectory <- "/rdss/ssa2724/fsmresfiles/People/collaborations/SubiaAhmed_HanxiaoSun/scRNAseq/Dec2023_MP2-E7_scRNAseq/Seurat/SA2_078S/"

# Initialize the Seurat object with the raw (non-normalized data).
SA2_078S <- CreateSeuratObject(counts = SA2_078S.data, project = "SA2_078", min.cells = 3, min.features = 200)

SA2_078S

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
SA2_078S[["percent.mt"]] <- PercentageFeatureSet(SA2_078S, pattern = "^MT-")

# Visualize QC metrics as a violin plot
#not necessary
ggsave(paste0(plotDirectory, "QC_vln_plot.png"), plot, width = 10, height = 8, units = "in")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(SA2_078S, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SA2_078S, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

SA2_078S <- subset(SA2_078S, subset = percent.mt < 5)

#SA2_078S <- subset(SA2_078S, subset = nFeature_RNA > 3000 & nFeature_RNA < 10000 & percent.mt < 5)

SA2_078S <- NormalizeData(SA2_078S, normalization.method = "LogNormalize", scale.factor = 10000)

SA2_078S <- NormalizeData(SA2_078S)

SA2_078S <- FindVariableFeatures(SA2_078S, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SA2_078S), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SA2_078S)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# save instead!!

all.genes <- rownames(SA2_078S)
SA2_078S <- ScaleData(SA2_078S, features = all.genes)

SA2_078S <- RunPCA(SA2_078S, features = VariableFeatures(object = SA2_078S))

# Examine and visualize PCA results a few different ways
print(SA2_078S[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(SA2_078S, dims = 1:2, reduction = "pca")

DimPlot(SA2_078S, reduction = "pca") + NoLegend()

ElbowPlot(SA2_078S, ndims = 30)

SA2_078S <- FindNeighbors(SA2_078S, dims = 1:25)
SA2_078S <- FindClusters(SA2_078S, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(SA2_078S), 5)

SA2_078S <- RunUMAP(SA2_078S, dims = 1:25)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(SA2_078S, reduction = "umap")

saveRDS(SA2_078S, file = "SA2_078S.rds")

# find all markers of cluster 4
cluster4.markers <- FindMarkers(SA2_078S, ident.1 = 4)
head(cluster4.markers, n = 5)

cluster2.markers <- FindMarkers(SA2_078S, ident.1 = 2)
head(cluster2.markers, n = 5)

VlnPlot(SA2_078S, features = c("PSG5", "BAHCC1"))
