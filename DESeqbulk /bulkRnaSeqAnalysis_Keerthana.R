#Code to analyse bulk RNAseq output from STAR - created by Keerthana on 2024-06-20 and last modified on 2024-06-20 at 2153
library(Seurat)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(ggplot2)


pathToFiles = "/projects/b1042/GoyalLab/aleona/STARAlignmentPipeline/Keerthana_outputnew/"
pathToFiles = "/projects/b1042/GoyalLab/aleona/DESEQ/rawData/bulkRNASeq/"
extData = "/projects/b1042/GoyalLab/aleona/DESEQ/extData/Keerthana_outputnew/"
extData = "/projects/b1042/GoyalLab/aleona/DESEQ/extData/rawdata1/"
files <- dir(path = pathToFiles,
             pattern = "*ReadsPerGene.out.tab", full.names = T) 

counttablefull <- files %>%
  map(read_tsv,  skip = 4, col_names = FALSE ) %>%
  purrr::reduce(cbind) 

datasets <-
  files %>%
  stringr::str_replace(paste0(pathToFiles,"/"), "") %>% # replace data/ at the beginning of filename
  stringr::str_replace("_ReadsPerGene.out.tab", "") # replace ReadsPerGene.out.tab end of name
datasets

columnnames <- c()
for (i in datasets) {
  columnnames <- c(columnnames,
                  paste0("gene_", i),
                  paste0(i, "_unstranded"),
                  paste0(i, "_forwardstrand"),
                  paste0(i, "_reversestrand")
  )
}


names(counttablefull) <- columnnames
rm(columnnames, datasets)

counttablefull <- counttablefull %>%
  mutate(ensembl_gene = "gene_ID:SA_MP2_E7_naive_1_S1_L001ReadsPerGene.out.tab") %>%
  dplyr::select(-starts_with("gene"))

counttablefull %>% head()


strandWiseCounts <- counttablefull %>%
                  dplyr::select(-ensembl_gene) %>%  # Use dplyr::select to ensure the right function is used
                  summarise(across(everything(), sum)) %>%
                  pivot_longer(cols = everything(), names_to = "library", values_to = "counts") %>%
                  separate(col = library, into = c("combined", "stranding"), sep = "_(?=[^_]+$)", extra = "merge") %>%
                  separate(col = combined, into = c("treatment", "replicate"), sep = "_", extra = "merge") %>%
                  pivot_wider(names_from = stranding, values_from = counts) %>%
                  mutate(
                    propF = forwardstrand / unstranded,  # Calculate the proportion of forwardstrand to unstranded
                    propR = reversestrand / unstranded   # Calculate the proportion of reversestrand to unstranded
                  )

write_csv(counttablefull,paste0(extData,"BulkRNAseq_countTableAll.csv"))
counttablefull <- read_csv(paste0(extData,"BulkRNAseq_countTableAll.csv"))

counttablefull <- counttablefull %>% 
  dplyr::select(-ends_with("_reversestrand"), -ends_with("forwardstrand")) 
names(counttablefull) <- str_replace(names(counttablefull), "_unstranded", "")

dim(counttablefull)

# Load the data using read.table from base R
# #This did not work because the ID is ENST for transcripts whereas earlier I had ENSG for genes. One gene can give many many transcripts - so not 1-1 mapping
# geneName <- read.table("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/bulkRNAseqData/genome/gencode.v46.metadata.HGNC",
#                        header = FALSE,  # Indicate no header in the file
#                        sep = "",        # Let R handle whitespace delimiters
#                        col.names = c("TranscriptID", "GeneName", "HGNC_ID"),  # Assign column names
#                        quote = "")      # Disable handling of quotes
# 
# library(biomaRt)
# listEnsembl("GRCh=38")

# geneCommon <- data.frame(geneName = intersect(check$GeneName, Features(naive2)))

geneFromGTF <- read_csv("/home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/bulkRNAseqData/genome/geneNamesFromGTF-2020.csv")
geneFromGTF$geneNameFinal <- sub("\\..*$", "", geneFromGTF$gene_id)

newJoinedData <- left_join(counttablefull, geneFromGTF, by = c("ensembl_gene" = "geneNameFinal"))
rownames(newJoinedData) <- newJoinedData$gene_id

bulkSeqData <- newJoinedData %>% select(-c("...1", "ensembl_gene", "gene_name", "SA2_028_MP2_E7_Gemcitabine_2_S13")) %>%
  as.data.frame() 
rownames(bulkSeqData) <- bulkSeqData$gene_id
bulkSeqData <- bulkSeqData %>% select(-c("gene_id"))

##
bulk_SeqData <- as.data.frame(counttablefull)
rownames(bulk_SeqData) <- bulk_SeqData$ensembl_gene
bulk_SeqData<- bulk_SeqData %>% dplyr::select(-c("ensembl_gene")) 
  
columnDataForDeDeq <- data.frame(sample = colnames(bulk_SeqData))
columnDataForDeDeq <- columnDataForDeDeq %>%
  mutate(condition = case_when(
    grepl("Gemcitabine", sample) ~ "Gemcitabine",
    grepl("Sotorasib", sample) ~ "Sotorasib",
    grepl("naive", sample) ~ "Naive",
    grepl("Cisplatin", sample) ~"Cisplatin"
  ))
rownames(columnDataForDeDeq) <- columnDataForDeDeq$sample
columnDataForDeDeq$type <- "paired-end"
# columnDataForDeDeq <- columnDataForDeDeq %>% select(-c("sample"))

# Assuming 'condition' is a character vector that needs to be converted to a factor
columnDataForDeDeq$condition <- factor(columnDataForDeDeq$condition,
                                       levels = c("Naive", "Gemcitabine", "Sotorasib", "Cisplatin"))

columnDataForDeDeq <- columnDataForDeDeq %>%
  group_by(condition) %>%  # Group by the 'condition' column
  mutate(replicate = row_number()) %>%  # Add a 'replicate' column with sequential numbering
  ungroup() %>%
  as.data.frame() # Remove the grouping

rownames(columnDataForDeDeq) <- columnDataForDeDeq$sample
columnDataForDeDeq <- columnDataForDeDeq %>% dplyr::select(-c("sample"))


# Now create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = bulk_SeqData,
                              colData = columnDataForDeDeq,
                              design = ~ condition)

smallestGroupSize <- 1
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$condition <- relevel(dds$condition, ref = "Naive")


dds <- DESeq(dds)

res_naive_vs_sotarasib2 <- results(dds, contrast=c("condition", "Naive", "Sotorasib"))
res_naive_vs_sotarasib2 <- data.frame(res_naive_vs_sotarasib)
res_naive_vs_sotarasib2 <- res_naive_vs_sotarasib %>% filter(padj < 0.05)

res_naive_vs_gemcitabine <- results(dds, contrast=c("condition", "Naive", "Gemcitabine"))
res_naive_vs_gemcitabine <- data.frame(res_naive_vs_gemcitabine)
res_naive_vs_gemcitabine <- res_naive_vs_gemcitabine %>% filter(padj < 0.05)

res_naive_vs_cisplatin <- results(dds, contrast=c("condition", "Naive", "Cisplatin"))
res_naive_vs_cisplatin <- data.frame(res_naive_vs_cisplatin)
res_naive_vs_cisplatin <- res_naive_vs_cisplatin %>% filter(padj < 0.05)

#Comparing Sotarasib to others
dds_sotarasib <- DESeqDataSetFromMatrix(countData = bulkSeqData,
                                               colData = columnDataForDeDeq,
                                               design = ~ condition)

smallestGroupSize <- 1
keep <- rowSums(counts(dds_sotarasib) >= 10) >= smallestGroupSize
dds_sotarasib <- dds_sotarasib[keep,]

dds_sotarasib$condition <- relevel(dds_sotarasib$condition, ref = "Sotorasib")

dds_sotarasib <- DESeq(dds_sotarasib)

res_sotarasib_vs_naive <- results(dds_sotarasib, contrast=c("condition", "Sotorasib", "Naive"))
res_sotarasib_vs_naive <- data.frame(res_sotarasib_vs_naive)
res_sotarasib_vs_naive <- res_sotarasib_vs_naive %>% filter(padj < 0.05)

res_sotarasib_vs_gemcitabine <- results(dds_sotarasib, contrast=c("condition", "Sotorasib", "Gemcitabine"))
res_sotarasib_vs_gemcitabine <- data.frame(res_sotarasib_vs_gemcitabine)
res_sotarasib_vs_gemcitabine <- res_sotarasib_vs_gemcitabine %>% filter(padj < 0.05)

res_sotarasib_vs_cisplatin <- results(dds_sotarasib, contrast=c("condition", "Sotorasib", "Cisplatin"))
res_sotarasib_vs_cisplatin <- data.frame(res_sotarasib_vs_cisplatin)
res_sotarasib_vs_cisplatin <- res_sotarasib_vs_cisplatin %>% filter(padj < 0.05)

#Comparing Gemcitabine to others
dds_gemcitabine <- DESeqDataSetFromMatrix(countData = bulkSeqData,
                                        colData = columnDataForDeDeq,
                                        design = ~ condition)

smallestGroupSize <- 1
keep <- rowSums(counts(dds_gemcitabine) >= 10) >= smallestGroupSize
dds_gemcitabine <- dds_gemcitabine[keep,]

dds_gemcitabine$condition <- relevel(dds_gemcitabine$condition, ref = "Gemcitabine")
dds_gemcitabine <- DESeq(dds_gemcitabine)

res_gemcitabine_vs_naive <- results(dds_gemcitabine, contrast=c("condition", "Gemcitabine", "Naive"))
res_gemcitabine_vs_naive <- data.frame(res_gemcitabine_vs_naive)
res_gemcitabine_vs_naive <- res_gemcitabine_vs_naive %>% filter(padj < 0.05)

res_gemcitabine_vs_sotarasib <- results(dds_gemcitabine, contrast=c("condition", "Gemcitabine", "Sotorasib"))
res_gemcitabine_vs_sotarasib <- data.frame(res_gemcitabine_vs_sotarasib)
res_gemcitabine_vs_sotarasib <- res_gemcitabine_vs_sotarasib %>% filter(padj < 0.05)

res_gemcitabine_vs_cisplatin <- results(dds_gemcitabine, contrast=c("condition", "Gemcitabine", "Cisplatin"))
res_gemcitabine_vs_cisplatin <- data.frame(res_gemcitabine_vs_cisplatin)
res_gemcitabine_vs_cisplatin <- res_gemcitabine_vs_cisplatin %>% filter(padj < 0.05)

# Function to add gene IDs as a column and prepare for merging
prepare_results <- function(df, comparison_name) {
  df$ensembl_id <- rownames(df)  # Add a new column for gene IDs from rownames
  df$comparison <- comparison_name  # Add a column for the comparison type
  rownames(df) <- NULL  # Optionally remove the rownames
  df <- merge(df, geneFromGTF, by.x = "ensembl_id", by.y = "gene_id", all.x = TRUE)
  return(df)
}

# Apply the function to each result set
res_naive_vs_sotarasib_final <- prepare_results(res_naive_vs_sotarasib, "Naive_vs_Sotorasib")
res_naive_vs_gemcitabine_final <- prepare_results(res_naive_vs_gemcitabine, "Naive_vs_Gemcitabine")
res_naive_vs_cisplatin_final <- prepare_results(res_naive_vs_cisplatin, "Naive_vs_Cisplatin")
res_sotarasib_vs_naive_final <- prepare_results(res_sotarasib_vs_naive, "Sotorasib_vs_Naive")
res_sotarasib_vs_gemcitabine_final <- prepare_results(res_sotarasib_vs_gemcitabine, "Sotorasib_vs_Gemcitabine")
res_sotarasib_vs_cisplatin_final <- prepare_results(res_sotarasib_vs_cisplatin, "Sotorasib_vs_Cisplatin")
res_gemcitabine_vs_naive_final <- prepare_results(res_gemcitabine_vs_naive, "Gemcitabine_vs_Naive")
res_gemcitabine_vs_sotarasib_final <- prepare_results(res_gemcitabine_vs_sotarasib, "Gemcitabine_vs_Sotorasib")
res_gemcitabine_vs_cisplatin_final <- prepare_results(res_gemcitabine_vs_cisplatin, "Gemcitabine_vs_Cisplatin")

# Combine all results into one data frame
all_results <- do.call(rbind, list(
  res_naive_vs_sotarasib_final,
  res_naive_vs_gemcitabine_final,
  res_naive_vs_cisplatin_final,
  res_sotarasib_vs_gemcitabine_final,
  res_sotarasib_vs_cisplatin_final,
  res_gemcitabine_vs_cisplatin_final
))

write_csv(all_results, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/deSeqBulkRnaSeq_withoutOneGem.csv")

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
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
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/bulkPCA_VST_withoutOneGem.png", pcaPlotVST)
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/bulkPCA_VST_withoutOneGem.svg", pcaPlotVST)

write_csv(pcaData_vst, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plotData/bulkPCA_VST__withoutOneGem.csv")

# pcaPlotNTD <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   theme_minimal()+
#   theme(legend.position = "bottom",  # Adjust this to put legend where you prefer
#         plot.title = element_text(hjust = 0.5),
#         plot.background = element_blank(),
#         panel.background = element_rect(fill = "white", colour = "black"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#         panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', colour = "gray"),
#         panel.grid.minor = element_blank(),
#         aspect.ratio = 0.6)+
#   coord_fixed()
# pcaPlotNTD
# ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/bulkPCA_NTD.png", pcaPlotNTD)
# write_csv(pcaData, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plotData/bulkPCA_NTD.csv")

BiocManager::install("vsn")
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

# Group by gene_name and comparison, then slice the first entry
all_results_heatmap <- all_results %>%
  distinct(gene_name, comparison, .keep_all = TRUE)


heatmap_data <- all_results_heatmap %>%
  select(gene_name, comparison, log2FoldChange) %>%
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

