#Code to analyze Count Input Matrix using DESeq2 - created by Aurelia on 2024-08-09 and is based on Keerthana's code 
# last modified on 03 Sep 2024
library(Seurat)
library(dplyr)
library(dbplyr)
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(readr)


#Main ---------------------------------
extData <- "/projects/b1042/GoyalLab/aleona/DESeq/extData/bulkRNASeq2/"
Scripts<- "/projects/b1042/GoyalLab/aleona/DESEQ/Script/"

source(paste0(Scripts, "deseqbulkfunction.R"))

# extData2 <- "/projects/b1042/GoyalLab/aleona/DESeq/extData/rawData2/"

countmatrix <- read_csv(paste0(extData, "countmatrix.csv"))
#ensembl_gene
countmatrixdf <- column_to_rownames(as.data.frame(countmatrix),'ensembl_gene')
countmatrixdf<- dplyr::select(countmatrixdf,-c("Cisplatin_1_S7", "Cisplatin_2_S8", "Gemcitabine_3_S14", "Gemcitabine_1_S12"))
# cts <- as.matrix.data.frame(countmatrixdf)

#Make an annotation matrix 
# First, let's create the annotation matrix
sample_names <- colnames(countmatrixdf)[colnames(countmatrixdf) != "ensembl_gene"]
annotation_data <- create_annotation_matrix(sample_names, patterns = "^(.+?)_.*")

# View the first few rows of the annotation matrix
head(annotation_data)

# annotate_data <- column_to_rownames(as.data.frame(annotate_data),'sample') 
# annotate_data <- annotate_data %>%
#   filter(!(condition == "Gemcitabine" & replicate == 2))
# annotation_data$condition <- factor(annotation_data$condition, levels = c("naive", "Gemcitabine", "Sotorasib", "Cisplatin"))
annotation_data$condition <- factor(annotation_data$condition, levels = c("naive", "Sotorasib"))
write_csv(annotation_data, paste0(extData,"annotationmatrix.csv"))

# annotate_data <- annotate_data %>%
#   group_by(condition) %>%  # Group by the 'condition' column
#   # mutate(replicate = row_number()) %>%  # Add a 'replicate' column with sequential numbering
#   ungroup() %>%
#   as.data.frame() # Remove the grouping

dds <- DESeqDataSetFromMatrix(countData = countmatrixdf,
                              colData = annotation_data,
                              design = ~ condition)


#Pre-filtering process-----------------------
smallestGroupSize <- 1
count_threshold <- 6

# Get the count matrix once
# Use base R's rowSums, which is still quite efficient
keep <- rowSums(counts(dds) >= count_threshold) >= smallestGroupSize

# Use logical indexing
dds <- dds[keep,]

##Differential Expression Analysis--------------------------------
# dds_naive <- dds
# dds_naive$condition <- relevel(x = dds_naive$condition, ref = "naive")
# 
# dds_naive <- DESeq(dds_naive)
# res_naive_vs_sotorasib <- results(dds_naive, contrast=c("condition","naive","Sotorasib"))
# res_naive_vs_sotarasib <- data.frame(res_naive_vs_sotorasib)
# res_naive_vs_sotarasib <- res_naive_vs_sotarasib %>% filter(padj < 0.05)

#Naive samples
contrast_info <- c("condition", "Sotorasib", "naive")
padj_threshold <- 0.05

# Relevel the condition factor to indicate that the reference is naive
dds$condition <- relevel(dds$condition, ref = "naive")

# Perform DESeq analysis
dds <- DESeq(dds)

# Get results for the specified comparison
full_results <- results(dds, contrast = contrast_info)
saveRDS(full_results, paste0(extData, "DESeqResults0.1_nvS.rds"))


#Adding gene_id to ENSEMBL only data and combine--------------------------------
ensembl.ids <- countmatrix['ensembl_gene']

ensembl_98 <- useEnsembl(biomart = 'genes', 
                          dataset = 'hsapiens_gene_ensembl',
                          version = 98)

gene_id <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
      filters = "ensembl_gene_id",
      values = ensembl.ids$ensembl_gene,
      mart = ensembl_98)

# Assuming you have your filtered results and gene_id dataframe ready
filtered_results_list <- list(full_results)
comparisons <- c("Naive_vs_Sotorasib")

all_results <- combine_deseq_results(filtered_results_list, gene_id, comparisons)


##Match Protein Classes ------------------------------------------
#Read the Protein atlast database 
proteinatlas <- read_tsv("/projects/b1042/GoyalLab/aleona/DESeq/Protein_data/proteinatlas.tsv", 
                         col_names = TRUE, 
                         show_col_types = FALSE)

# Assuming proteinatlas is already loaded

full_results <- match_and_classify_proteins(all_results, "ensembl_gene", proteinatlas)

# Ordered by abs lfc value
ordered <- as.data.frame(full_results)

# Convert to absolute values
ordered$log2FoldChange <- abs(ordered$log2FoldChange)

# Order the data frame
ordered <- ordered[order(ordered$log2FoldChange, decreasing = TRUE), ]

#Only the Log2FoldChange above 2 
saveabove2 <- ordered %>% filter(log2FoldChange > 2 ) %>% filter(padj < 0.05 )

#Saving results 
# write_csv(saveabove2, paste0(extData2, "keeptrackSep32024.csv"))
write_csv(all_results, paste0(extData, "DESeqwithgeneid_NaivevsSotfull2.csv"))
write_csv(saveabove2, paste0(extData, "DESeq_NvSabsLFCresult.csv"))
saveRDS(dds, paste0(extData, "DESEQ.rds"))

#Filter the data------------------------------------------
ordered1 <- as.data.frame(full_results)

#Split between overexpression and underexpression 
over <- ordered1 %>% filter(log2FoldChange > 2 ) %>% filter(padj < 0.05 ) 
under <- ordered1 %>% filter(log2FoldChange < -2 ) %>% filter(padj < 0.05 )

#Get the gene count data or rawdata
genecountdata <- countmatrixdf
genecountdata <- rownames_to_column(genecountdata, "ensembl_gene")

#Add the rawcount into the dataset 
#For left_join to work, you need to get the data 
over_all <- left_join(over, genecountdata, by = "ensembl_gene")
under_all <- left_join(under, genecountdata, by = "ensembl_gene")

write_csv(over_all, paste0(extData, "overexpress.csv"))
write_csv(under_all, paste0(extData, "underexpress.csv"))

#Filter just surfaceproteins 
over_sprotein <- over_all %>% filter(`Protein Class`=="Surface Proteins")
under_sprotein <- under_all %>% filter(`Protein Class`=="Surface Proteins")

write_csv(over_sprotein, paste0(extData, "overexpress_surfaceprot.csv"))
write_csv(under_sprotein, paste0(extData, "underexpress_surfaceprot.csv"))


# saveabove2 <- read_csv("/projects/b1042/GoyalLab/aleona/DESeq/extData/DESeq_NvSabsLFCresult.csv")
# 
# #Making into dataframe 
# dataset <- as.data.frame(saveabove2)
# 
# #Making into only surface proteins 
# surfaceproteins <- dataset %>% filter(`Protein Class`=="Surface Proteins")
# 
# rawdata <-countmatrixdf 
# 
# filtered_df <- rawdata %>%
#   rownames_to_column("ensembl_gene") %>%
#   filter(ensembl_gene %in% surfaceproteins$ensembl_gene)



# ##Introduction to DGE------------------------------------- 
# #https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html
# res_tableOE <- lfcShrink(dds, coef =3, res= full_results)
# summary(res_tableOE)
# 
# ### Set thresholds
# padj.cutoff <- 0.05
# lfc.cutoff <- 0.58
# 
# res_tableOE_tb <- res_tableOE %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble()
# 
# sigOE <- res_tableOE_tb %>%
#   filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

# 
# 
# ##Independent hypothesis weighting------------------------------------- 
# library("IHW")
# dds_naive <- dds
# dds_naive$condition <- relevel(x = dds_naive$condition, ref = "naive")
# 
# dds_naive <- DESeq(dds_naive)
# res <- results(dds_naive, contrast=c("condition","naive","Sotorasib"))
# res_naive_vs_sotarasib <- data.frame(res) %>% filter(padj < 0.05)
# resnvSIHW <- results(dds_naive, contrast=c("condition","naive","Sotorasib"), filterFun=ihw)
# summary(resnvSIHW)
# sum(resnvSIHW$padj < 0.05, na.rm=TRUE)
# metadata(resnvSIHW)$ihwResult
# 
# 
# 
# #Testing between 0.05 and 0.1-------------------------------------------
# dds_test <- DESeq(dds)
# res_05 <- results(dds_test, alpha=0.05)
# res_10 <- results(dds_test, alpha=0.1)
# 
# # Compare log2FoldChange values
# summary(res_05$log2FoldChange - res_10$log2FoldChange)
# 
# # Compare lfcSE values
# summary(res_05$lfcSE - res_10$lfcSE)
# 
# # Compare p-values
# summary(res_05$pvalue - res_10$pvalue)
# 
# # Compare adjusted p-values
# summary(res_05$padj - res_10$padj)
# 
# # Compare adjusted p-values
# summary(res_05$padj - res_10$padj)
