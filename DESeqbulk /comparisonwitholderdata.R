# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)  # For adding text labels

# Read the two CSV files
extData <- "/projects/b1042/GoyalLab/aleona/DESEQ/extData/bulkRNASeq2/"
data1 <- read_csv("/projects/b1042/GoyalLab/aleona/DESEQ/extData/bulkRNASeq2/DESeqwithgeneid_NaivevsSotfull2.csv")
# data2 <- read_csv("/projects/b1042/GoyalLab/aleona/DESEQ/KeerthanaDESEQ/deSeqBulkRnaSeq_withoutOneGem.csv")
data2 <- read_csv("/projects/b1042/GoyalLab/aleona/DESEQ/extData/rawData2/DESeqwithgeneid_NaivevsSotfull2.csv")
data1f <- filter(data1, padj<=0.05)
data2f <- filter(data2, comparison=="Naive_vs_Sotorasib")
data2f <- filter(data2f, padj<=0.05)
# Ensure both datasets have a common identifier column (e.g., gene_id)
# If not already present, you might need to create one

# # Merge the two datasets
# merged_data <- full_join(
#   data1f %>% 
#     dplyr::select(ensembl_gene, log2FoldChange, padj) %>% 
#     dplyr::rename(log2FC1 = log2FoldChange, padj1 = padj),
#   data2f %>% 
#     dplyr::select(geneNameFinal, log2FoldChange, padj) %>% 
#     dplyr::rename(log2FC2 = log2FoldChange, padj2 = padj, ensembl_gene = geneNameFinal),
#   by = "ensembl_gene"
# )

merged_data <- full_join(
  data1f %>% 
    dplyr::select(ensembl_gene, log2FoldChange, padj) %>% 
    dplyr::rename(log2FC1 = log2FoldChange, padj1 = padj),
  data2f %>% 
    dplyr::select(ensembl_gene, log2FoldChange, padj) %>% 
    dplyr::rename(log2FC2 = log2FoldChange, padj2 = padj),
  by = "ensembl_gene"
)

# Create a scatter plot comparing log2FoldChanges
comparison_plot <- ggplot(merged_data, aes(x = log2FC1, y = log2FC2)) +
  geom_point(aes(color = (padj1 < 0.05 & padj2 < 0.05)), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("grey", "blue")) +
  labs(
    title = "Comparison of Log2 Fold Changes",
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant in both"
  ) +
  theme_minimal() +
  coord_fixed()  # Ensure x and y scales are the same

# Add labels for some interesting points (e.g., highly discordant genes)
interesting_genes <- merged_data %>%
  mutate(diff = abs(log2FC1 - log2FC2)) %>%
  filter(padj1 < 0.05 & padj2 < 0.05) %>%
  top_n(5, diff)

comparison_plot <- comparison_plot +
  geom_text_repel(
    data = interesting_genes,
    aes(label = ensembl_gene),
    box.padding = 0.5,
    max.overlaps = Inf
  )

# Display the plot
print(comparison_plot)

# Save the plot
ggsave("log2FC_comparison_plotfin.png", comparison_plot, width = 10, height = 8, dpi = 300)

# Calculate correlation
correlation <- cor(merged_data$log2FC1, merged_data$log2FC2, use = "complete.obs")
print(paste("Correlation between log2FoldChanges:", round(correlation, 3)))

# Summarize discordant genes (significant in one dataset but not the other)
discordant_genes <- merged_data %>%
  filter((padj1 < 0.05 & padj2 >= 0.05) | (padj1 >= 0.05 & padj2 < 0.05))

print(paste("Number of discordant genes:", nrow(discordant_genes)))

# Output discordant genes to a CSV file
write_csv(discordant_genes, paste0(extData,"discordant_genes.csv"))

genes_only_in_data2 <- anti_join(
  data2f %>% dplyr::select(geneNameFinal ) %>% dplyr::rename(ensembl_gene = geneNameFinal),
  data1f %>% dplyr::select(ensembl_gene),
  by = "ensembl_gene"
)

genes_only_in_data2 <- anti_join(
  data2f %>% dplyr::select(ensembl_gene),
  data1f %>% dplyr::select(ensembl_gene),
  by = "ensembl_gene"
)


# Now, let's create a new dataframe with the genes only in data1f 
# and their corresponding log2FoldChange values from both datasets
result <- genes_only_in_data2 %>%
  left_join(
    data2f %>% 
      dplyr::select(ensembl_gene, log2FC2 = log2FoldChange, padj2=padj, pvalue2=pvalue),
    by = "ensembl_gene"
  )  %>%
  left_join(
    data1 %>% 
      dplyr::select(ensembl_gene, log2FC1 = log2FoldChange, padj1=padj, pvalue1=pvalue),
    by = "ensembl_gene"
  ) 

result <- genes_only_in_data2 %>%
  left_join(genes_only_in_data2_with_fc, by = "ensembl_gene")

genes_only_in_data1 <- anti_join(
  data1f %>% dplyr::select(ensembl_gene),
  data2f %>% dplyr::select(ensembl_gene),
  by = "ensembl_gene"
)

result_genedata1 <- genes_only_in_data1 %>%
  left_join(
    data1f %>% 
      dplyr::select(ensembl_gene, log2FC1 = log2FoldChange, padj1=padj, pvalue1=pvalue),
    by = "ensembl_gene"
  )  %>%
  left_join(
    data2 %>% 
      dplyr::select(ensembl_gene, log2FC2 = log2FoldChange, padj2=padj, pvalue2=pvalue),
    by = "ensembl_gene"
  ) 

##See if there are duplicates ------------------------
sum(duplicated(data1f$ensembl_gene))
sum(duplicated(data2f$geneNameFinal))

duplicates_data2 <- data2f[duplicated(data2f$geneNameFinal) | duplicated(data2f$geneNameFinal, fromLast = TRUE), ]
print(duplicates_data2)

##Making a Venn Diagram ---------------
if (!requireNamespace("VennDiagram", quietly = TRUE))
  install.packages("VennDiagram")
library(VennDiagram)
library(grid)

genes_data1 <- data1$ensembl_gene
genes_data2 <- data2$ensembl_gene

# Clear all graphics devices
graphics.off()
while (!is.null(dev.list()))  
  dev.off()
dev.new(width=10, height=10)


# Create your Venn diagram object
venn.plot <- draw.pairwise.venn(
  area1 = length(unique(genes_data1)),
  area2 = length(unique(genes_data2)),
  cross.area = length(intersect(genes_data1, genes_data2)),
  category = c("", ""),  # We'll add custom labels later
  fill = c("lightblue", "lightgreen"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(0, 0),  # Center the diagrams
  cat.dist = 0.05,
  cat.just = list(c(-1, -1), c(1, 1)),
  ext.pos = 30,
  ext.dist = -0.05,
  ext.length = 0.85,
  ext.line.lwd = 2,
  ext.line.lty = "dashed"
)

# Function to create a legend
create_legend <- function() {
  legend_plot <- gTree(children = gList(
    rectGrob(x = c(0.1, 0.1), y = c(0.7, 0.3), width = 0.05, height = 0.05, 
             gp = gpar(fill = c("lightblue", "lightgreen"), col = "black")),
    textGrob(c("All Data", "Naive and Sotorasib"), x = 0.2, y = c(0.7, 0.3), just = "left",
             gp = gpar(fontsize = 14))
  ))
  return(legend_plot)
}

# Create a new plot window
grid.newpage()

# Create a viewport layout
pushViewport(viewport(layout = grid.layout(1, 2, widths = c(4, 1))))

# Draw Venn diagram
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(venn.plot)
popViewport()

# Draw legend
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(create_legend())
popViewport()

dev.off()
genes_in_both <- intersect(genes_data1, genes_data2)
genes_only_in_data1 <- setdiff(genes_data1, genes_data2)
genes_only_in_data2 <- setdiff(genes_data2, genes_data1)

cat("Unique genes in Data1:", length(unique(genes_data1)), "\n")
cat("Unique genes in Data2:", length(unique(genes_data2)), "\n")
cat("Genes in both datasets:", length(genes_in_both), "\n")
cat("Genes only in Data1:", length(genes_only_in_data1), "\n")
cat("Genes only in Data2:", length(genes_only_in_data2), "\n")




