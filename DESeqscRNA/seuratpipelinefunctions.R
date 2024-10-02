#This is the functions used in all the process of Seurat pipeline for scRNA sequencing 
# library(Seurat)
# library(harmony)
# library(dplyr)
# library(ggplot2)
# library(cowplot)
# library(viridis) 
# library(RColorBrewer)
# library(reticulate)
# library(readr)
# library(ggrepel)
# library(patchwork)
# library(R.utils)  # For gunzip function

# Load required libraries
required_packages <- c(
  "Seurat", "harmony", "dplyr", "ggplot2", "cowplot", "viridis", "RColorBrewer", "reticulate", 
  "readr", "ggrepel", "patchwork", "R.utils"
  )

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)


#' Create Directory if It Doesn't Exist
#'
#' This function creates a directory at the specified path if it doesn't already exist.
#'
#' @param dir_path A character string specifying the path of the directory to be created.
#'
#' @return No return value. The function creates the directory and prints a message if successful.
#'
#' @examples
#' create_dir_if_not_exists("path/to/new/directory")
create_dir_if_not_exists <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    print(paste("Created directory:", dir_path))
  }
}

# Function to read and initialize Seurat objects from 10X data
# This follows the Seurat pipeline 
initSeuratObject <- function(directory, projectName) {
  #Read in the 10x data
  counts <- Read10X(directory)
  # Initialize the Seurat object with the raw (non-normalized data).
  seuratObject <- CreateSeuratObject(counts, min.cells = 3, min.features = 200, project = projectName)
  # QC metrics of percent mitochondria is stored here 
  seuratObject[["percent.mito"]] <- PercentageFeatureSet(seuratObject, pattern = "^MT-")
  
  return(seuratObject)
}


# Function to create and save QC plots
createAndSaveQCPlots <- function(seuratObject, path, prefix, percentageMitoCutoff, minimumCounts, minimumFeatures, maximumCounts = NA, maximumFeatures = NA) {
  # Scatter plot of RNA counts vs. mitochondrial percentage
  plot1 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "percent.mito") + 
    geom_vline(xintercept = minimumCounts, linetype = "dotted", color = "red") + 
    geom_hline(yintercept = percentageMitoCutoff, linetype = "dashed", color = "blue")
  if (!is.na(maximumCounts)) {
    plot1 <- plot1 + geom_vline(xintercept = maximumCounts, linetype = "dashed", color = "red")
  }
  
  # Scatter plot of RNA features vs. mitochondrial percentage
  plot2 <- FeatureScatter(seuratObject, feature1 = "nFeature_RNA", feature2 = "percent.mito") + 
    geom_vline(xintercept = minimumFeatures, linetype = "dotted", color = "red") + 
    geom_hline(yintercept = percentageMitoCutoff, linetype = "dashed", color = "blue")
  if (!is.na(maximumFeatures)) {
    plot2 <- plot2 + geom_vline(xintercept = maximumFeatures, linetype = "dashed", color = "red")
  }
  
  # Scatter plot of RNA counts vs. RNA features
  plot3 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(xintercept = minimumCounts, linetype = "dotted", color = "red") + 
    geom_hline(yintercept = minimumFeatures, linetype = "dotted", color = "red")
  if (!is.na(maximumCounts) && !is.na(maximumFeatures)) {
    plot3 <- plot3 + geom_vline(xintercept = maximumCounts, linetype = "dashed", color = "red") + 
      geom_hline(yintercept = maximumFeatures, linetype = "dashed", color = "red")
  }
  
  
  # Create individual violin plots
  plot_mito <- create_violin(seuratObject, "percent.mito", max_threshold = percentageMitoCutoff)
  plot_count <- create_violin(seuratObject, "nCount_RNA", min_threshold= minimumCounts)
  plot_feature <- create_violin(seuratObject, "nFeature_RNA", min_threshold= minimumFeatures, max_threshold = maximumFeatures)
  
  # Combine plots
  combined_plot <- plot_mito | plot_count | plot_feature
  
  # Add titles
  combined_plot <- combined_plot +
    plot_annotation(
      title = "Quality Control Metrics",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
    )
  
  # Save plots to specified path
  ggsave(paste0(path, prefix, "CountVsPercentMito.png"), plot = plot1)
  ggsave(paste0(path, prefix, "FeatureVsPercentMito.png"), plot = plot2)
  ggsave(paste0(path, prefix, "CountVsFeature.png"), plot = plot3)
  ggsave(paste0(path, prefix, "CombinedQCPlot.png"), plot = combined_plot, 
         width = 24, height = 8, dpi = 300, units = "in")
}


# Function to create a single violin plot with a threshold line
create_violin <- function(seuratObject, feature, min_threshold = NULL, max_threshold = NULL) {
  p <- VlnPlot(seuratObject, features = feature, pt.size = 0) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(5, 5, 5, 5))
  
  if (!is.null(min_threshold)) {
    if (feature == "percent.mito") {
      p <- p + geom_hline(yintercept = min_threshold, linetype = "dashed", color = "red")
    } else {
      p <- p + geom_hline(yintercept = min_threshold, linetype = "dashed", color = "red")
    }
  }
  
  if (!is.null(max_threshold)) {
    if (feature == "percent.mito") {
      p <- p + geom_hline(yintercept = max_threshold, linetype = "dashed", color = "blue")
    } else {
      p <- p + geom_hline(yintercept = max_threshold, linetype = "dashed", color = "blue")
    }
  }
  
  return(p)
}


#Function to plot QC plots for all samples in one go
plotAll <- function(params, qcPlotFilePath) {
  lapply(params, function(p) {
    createAndSaveQCPlots(p$seuratObject, qcPlotFilePath, p$prefix, 
                         p$percentMitoCutoff, p$minCounts, p$minFeatures,
                         p$maxCounts, p$maxFeatures)
  })
}

# Function to apply filtering for a specific dataset
applyFilter <- function(seuratObject, params) {
  minCounts = params$minCounts
  minGenes = params$minFeatures
  maxMitoPercentage = params$percentMitoCutoff
  maxCounts = params$maxCounts
  maxGenes = params$maxFeatures
  
  seuratObject <- subset(seuratObject, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes &
                           percent.mito < maxMitoPercentage & nCount_RNA > minCounts & nCount_RNA < maxCounts)
  
  return(seuratObject)
}

# Function to perform integration and normalization
integrate_and_normalize <- function(object_list, n_pcs = 50, path_to_folder, reduction_name="harmony") {
  # Merge all objects
  combined <- merge(object_list[[1]], 
                    y = object_list[-1], 
                    add.cell.ids = names(object_list), 
                    project = "Merged",
                    merge.data = TRUE)
  
  # Preprocess the combined object ( Normalize, feature selection, scale and runpca)
  combined <- NormalizeData(combined)
  combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
  combined <- ScaleData(combined)
  combined <- RunPCA(combined, npcs = n_pcs)
  
  # Create and save elbow plot
  elbow_plot <- ElbowPlot(combined, ndims = n_pcs)
  ggsave(paste0(path_to_folder, "elbow_plot.png"), plot = elbow_plot, width = 10, height = 8, dpi = 300)
  
  # Run integration
  integrated <- IntegrateLayers(
    object = combined, 
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = reduction_name,
    features = Features(combined)
  )
  
  integrated <- JoinLayers(integrated)
  
  return(integrated)
}


# Function to perform umap and plot umap 
runumap <- function(object, n_pcs = 50, path_to_folder, reduction_name="harmony", custom_name = NULL) {
  # Get the name of the object
  object_name <- deparse(substitute(object))
  
  # Create a timestamp
  timestamp <- format(Sys.time(), "%Y%m%d")
  
  # Construct the filename
  if (!is.null(custom_name)) {
    # Use the custom name if provided
    filename <- paste0(custom_name, "_", timestamp, "_umap.png")
  } else {
    # Use the default naming convention if no custom name is provided
    filename <- paste0(object_name, "_", timestamp,"_umap.png")
  }
  
  # Run UMAP and clustering on integrated object
  object <- FindNeighbors(object, reduction = reduction_name, dims = 1:n_pcs)
  object <- FindClusters(object, resolution = 0.5)
  object <- RunUMAP(object, reduction = reduction_name, dims = 1:n_pcs)
  
  
  # Create and save UMAP plot
  umap_plot <- DimPlot(object = object, reduction = "umap", pt.size = .1, group.by = "orig.ident")
  ggsave(paste0(path_to_folder, filename), plot = umap_plot, width = 16, height = 12, dpi = 300)
  
  return(object)
}

runfeatureplot <- function(object, list_features, path_to_folder, custom_name = NULL,
                           cols = c("lightgrey" ,"blue") )
  {
  # Get the name of the object
  object_name <- deparse(substitute(object))
  
  # Create a timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  
  # Construct the filename
  if (!is.null(custom_name)) {
    # Use the custom name if provided
    filename <- paste0(custom_name, "_", timestamp, "_featureplot.png")
  } else {
    # Use the default naming convention if no custom name is provided
    filename <- paste0(object_name, "_", timestamp, "_featureplot.png")
  }
  
  # Ensure the path ends with a slash
  if (!endsWith(path_to_folder, "/")) {
    path_to_folder <- paste0(path_to_folder, "/")
  }
  
  
  # Create and save feature plots with customized parameters
  feature_plots <- FeaturePlot(object, 
                               features = list_features, 
                               ncol = 3, 
                               cols=cols)   # Set custom colors
  
  # Modify the color scale for when plots is all 0s
  for (i in seq_along(feature_plots)) {
    current_plot <- feature_plots[[i]]$data
    feature_name <- list_features[[i]]
    if ( all(current_plot[[feature_name]] == 0) ) {
      # If all values are 0, set all points to light gray
      feature_plots[[i]] <- FeaturePlot(object, 
                                        features = feature_name, 
                                        ncol = 3, 
                                        cols=c(cols[[1]],cols[[1]]))
    }
  }

  # Full path for saving the file
  full_path <- paste0(path_to_folder, filename)
  
  # Save the plot
  ggsave(full_path, plot = feature_plots, width = 16, height = 12, dpi = 300)
  
  # Print a message confirming the save location
  print(paste("Feature plot saved as:", full_path))
  
  return(feature_plots)
}
