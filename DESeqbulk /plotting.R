#Code to plot DESeq2 data - created by Aurelia on 2024-08-09 and is based on Keerthana's code 
# last modified on 
library(Seurat)
library(dplyr)
library(dbplyr)
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(readr)



#Functions ------------------------------
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





#Main -------------
extData <- "/projects/b1042/GoyalLab/aleona/DESEQ/extData/bulkRNAseq/"
plotfile <- "/projects/b1042/GoyalLab/aleona/DESEQ/plots/bulkRNAseq/"
create_dir_if_not_exists(plotfile)


dds <- readRDS(paste0(extData, "DESEQ.rds"))

#Extracting transformed values 
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)


##PCA Plot ---------------------------------------------
plotPCA(vsd, intgroup=c("condition", "replicate"))
ntd <- normTransform(dds)
pcaData <- plotPCA(ntd, intgroup=c("condition", "type"), returnData=TRUE)
write_csv(pcaData, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/bulkRNAseq/pcaCoordinates.csv")


vst_data <- vst(dds, blind = T)
pcaData_vst <- plotPCA(vst_data, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData_vst, "percentVar"))

pcaPlotVST <- ggplot(pcaData_vst, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()+
  theme(legend.position = "bottom",  # Adjust this to put legend where you prefer
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', colour = "gray"),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.6)+
  coord_fixed()
pcaPlotVST



##Heatmap ------------------------------------------
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$replicate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

##HEatmap of comparisons------------------------------------------
all_results <- read_csv(paste0(extData, "DESeqwithgeneid_NaivevsSot.csv"))

all_results_heatmap <- all_results %>%
  distinct(external_gene_name, comparison, .keep_all = TRUE)

heatmap_data <- all_results_heatmap %>%
  dplyr::select(external_gene_name, comparison, log2FoldChange) %>%
  spread(key = comparison, value = log2FoldChange)  # Adjust this if you want to visualize another metric


# Fill NAs if necessary
heatmap_data[is.na(heatmap_data)] <- 0

pheatmap(as.matrix(heatmap_data[,-1]),  # Exclude the gene column for the heatmap
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE,
         show_colnames = TRUE)

gene_list <- unique(all_results_heatmap$ensembl_id)
filtered_counts <- bulkSeqData
# Example of log transformation with a pseudocount
log_counts <- log2(filtered_counts + 1)

geneList <- data.frame(gene = unique(newJoinedData$gene_name))
write_csv(geneList, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/bulkRNAseqGeneList.csv")


##Plot MA function of normal dataset--------------------------------------
nvS_0.05 <- readRDS(paste0(extData, "DESeqResults0.05_nvS.rds"))
plotMA(nvS_0.05)


nvS_0.1 <- readRDS(paste0(extData, "DESeqResults0.1_nvS.rds"))
plotMA(nvS_0.1)



##Log fold change shrinkage for visualization and ranking------
dds_copy <- dds

# Relevel the condition factor
dds_copy$condition <- relevel(dds_copy$condition, ref = "naive")

# Perform DESeq analysis
dds_copy <- DESeq(dds_copy)

# Get results for the specified comparison
res <- results(dds_copy, contrast = c("condition",  "naive", "Sotorasib"), alpha = 0.05)

resultsNames(dds_copy)
resLFC <- lfcShrink(dds_copy, coef="condition_Sotorasib_vs_naive", type="apeglm")


##ViolinPlot---------------------------

# Read the CSV file
data <- read_csv(paste0(extData,"DESeqwithgeneid_NaivevsSotfull2.csv"))


# Subset the data where padj < 0.05
data_filtered <- data %>% filter(padj < 0.05003)

# Create the violin plot with the filtered data and subtle scatter dots
violin_plot <- ggplot(data_filtered, aes(x = comparison, y = log2FoldChange, fill = comparison)) +
  # Add the violin plot
  geom_violin(trim = FALSE, alpha = 0.7) +  # Increase alpha for more prominent violin
  
  # Add a boxplot inside the violin
  geom_boxplot(width = 0.1, fill = "white", color = "black", alpha = 0.7, outlier.shape = NA) +  # Remove boxplot outliers
  
  # Add jittered points (less prominent)
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) +  # Reduced width, size, and alpha for subtlety
  
  # Set the theme to classic (no gridlines, white background)
  theme_classic() +
  
  # Add labels to the plot
  labs(title = "Distribution of Log2 Fold Change by Comparison (padj < 0.05)",
       subtitle = paste("Number of genes:", nrow(data_filtered)),
       x = "Comparison",
       y = "Log2 Fold Change") +
  
  # Customize the theme
  theme(
    legend.position = "none",  # Remove the legend
    plot.title = element_text(hjust = 0.5),  # Center the title
    plot.subtitle = element_text(hjust = 0.5),  # Center the subtitle
    axis.text.x = element_text(angle = 45, hjust = 1)  # Angle x-axis labels for readability
  ) +
  
  # Flip the coordinates to make the plot horizontal
  coord_flip()


# Display the plot
print(violin_plot)

ggsave(paste0(plotfile,"log2FC_violinplot.png"), violin_plot, width = 10, height = 8, dpi = 300)






