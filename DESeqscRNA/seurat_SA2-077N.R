# *************************************************
# Seurat script. started by Subia Ahmed on February 28, 2024 
# **************************************************

library(dplyr)
library(Seurat)
library(patchwork)

# Load the datasets
SA2_077N.data <- Read10X(data.dir = "/projects/b1042/GoyalLab/Subia/20230122_scRNAseq/run_count_SA2_077/outs/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
SA2_077N <- CreateSeuratObject(counts = SA_077N.data, project = "SA2_077", min.cells = 3, min.features = 200)

SA2_077N

# Lets examine a few genes in the first thirty cells
# SA_077N.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
SA2_077N[["percent.mt"]] <- PercentageFeatureSet(SA2_077N, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(SA2_077N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(SA2_077N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SA2_077N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

SA2_077N <- subset(SA2_077N, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 10)

SA2_077N <- NormalizeData(SA2_077N, normalization.method = "LogNormalize", scale.factor = 10000)

SA2_077N <- NormalizeData(SA2_077N)

SA2_077N <- FindVariableFeatures(SA2_077N, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SA2_077N), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SA2_077N)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# save instead!!

all.genes <- rownames(SA2_077N)
SA2_077N <- ScaleData(SA2_077N, features = all.genes)

SA2_077N <- RunPCA(SA2_077N, features = VariableFeatures(object = SA2_077N))

# Examine and visualize PCA results a few different ways
print(SA2_077N[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(SA2_077N, dims = 1:2, reduction = "pca")

DimPlot(SA2_077N, reduction = "pca") + NoLegend()

ElbowPlot(SA2_077N)

SA2_077N <- FindNeighbors(SA2_077N, dims = 1:15)
SA2_077N <- FindClusters(SA2_077N, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(SA2_077N), 5)

SA2_077N <- RunUMAP(SA2_077N, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(SA2_077N, reduction = "umap")

saveRDS(SA2_077N, file = "SA2_077N.rds")
