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
extData <- "/projects/b1042/GoyalLab/aleona/bulk/extractedData/SA_HS_MP2E7/matrix/"
Scripts<- "/projects/b1042/GoyalLab/aleona/LabPipelines/DESeqbulk/"
output <- "/projects/b1042/GoyalLab/aleona/bulk/extractedData/SA_HS_MP2E7/DESeq/"

source(paste0(Scripts, "deseqbulkfunction.R"))

#Create directory 
create_dir_if_not_exists(output)

countmatrix <- read_csv(paste0(extData, "countmatrix.csv"))

#Look at what the column names are
print(colnames(countmatrix))


#  [1] "ensembl_gene"             "ID_HS_017_B1_sot_S84"     "ID_HS_017_B2_sot_S52"     "ID_HS_017_B3_sot_S59"    
# [5] "ID_HS_017_gem1_S51"       "ID_HS_017_gem2_S83"       "ID_HS_017_p2a2_S85"       "ID_HS_025_tram1_S60"     
# [9] "ID_HS_025_tram2_S69"      "ID_HS_025_tram3_S77"      "ID_HS_026_20g_20s1_S62"   "ID_HS_026_20g_20s2_S93"  
# [13] "ID_HS_026_20g_50s_S80"    "ID_HS_026_50g_20s1_S48"   "ID_HS_026_50g_20s2_S49"   "ID_HS_026_50g_50s1_S67"  
# [17] "ID_HS_026_50g_50s2_S76"   "ID_HS_027_12daysS1_S65"   "ID_HS_027_12daysS2_S86"   "ID_HS_027_12daysS3_S96"  
# [21] "ID_HS_027_48hCis1_S56"    "ID_HS_027_48hCis2_S61"    "ID_HS_027_48hCis3_S88"    "ID_HS_027_48hGem1_S66"   
# [25] "ID_HS_027_48hGem2_S73"    "ID_HS_027_48hGem3_S75"    "ID_HS_027_48hSot1_S64"    "ID_HS_027_48hSot2_S78"   
# [29] "ID_HS_027_48hSot3_S95"    "ID_HS_027_8daysG1_S71"    "ID_HS_027_8daysG2_S92"    "ID_HS_027_8daysG3_S54"   
# [33] "ID_HS_027_cis1_S58"       "ID_HS_027_cis2_S68"       "ID_HS_027_cis3_S90"       "ID_SA2-100_20g_20s_S79"  
# [37] "ID_SA2-100_20g_50s_S47"   "ID_SA2-100_50g_50s_S82"   "ID_SA2-101_50g_20s_S81"   "ID_SA2-102_A1_sotGem_S91"
# [41] "ID_SA2-102_A2_sotGem_S70" "ID_SA2-102_C1_sotR_S45"   "ID_SA2-102_C2_sotR_S57"   "ID_SA2-102_C3_sotR_S89"  
# [45] "ID_SA2-102_P_S53"         "ID_SA2-103_A2_sotGem_S63" "ID_SA2-103_P1_S74"        "ID_SA2-104_A1_gemSot_S94"
# [49] "ID_SA2-104_A2_gemSot_S72" "ID_SA2-104_B2_gem_S50"    "ID_SA2-104_C1_gemR_S55"   "ID_SA2-104_C2_gemR_S87"  
# [53] "ID_SA2-104_C3_gemR_S46"  

# Make ensembl_gene as the rownames for later processing 
countmatrixdf <- column_to_rownames(as.data.frame(countmatrix),'ensembl_gene')

# Select the columns you would like to have in your DESeq 
# if you want to delete the columns use -c("colname1", "colname2")
# If you want to include the columns use c("colname1", "colname2")
countmatrixdf<- dplyr::select(countmatrixdf,c("ID_HS_017_B1_sot_S84", "ID_HS_017_B2_sot_S52", "ID_HS_017_B3_sot_S59", "ID_HS_017_p2a2_S85","ID_SA2-102_P_S53", "ID_SA2-103_P1_S74"))
# countmatrixdf<- dplyr::select(countmatrixdf,-c("Cisplatin_1_S7", "Cisplatin_2_S8", "Gemcitabine_3_S14", "Gemcitabine_1_S12"))
# cts <- as.matrix.data.frame(countmatrixdf)

#Make an annotation matrix 
# First, let's create the annotation matrix
sample_names <- colnames(countmatrixdf)[colnames(countmatrixdf) != "ensembl_gene"]
# Specify the condition and the replicate 
condition_maps <- list(
  "ID_HS_017_B1_sot_S84" = list(condition = "Sotorasib", replicate = 1),
  "ID_HS_017_B2_sot_S52" = list(condition = "Sotorasib", replicate = 2),
  "ID_HS_017_B3_sot_S59" = list(condition = "Sotorasib", replicate = 3),
  "ID_HS_017_p2a2_S85" = list(condition = "Naive", replicate = 1),
  "ID_SA2-102_P_S53" = list(condition = "Naive", replicate = 2),
  "ID_SA2-103_P1_S74" = list(condition = "Naive", replicate = 3)
)
# Use the function to create annotation matrix 
annotation_data <- create_annotation_matrix(sample_names, condition_mapping=condition_maps)

# View the first few rows of the annotation matrix
head(annotation_data)

# Create a factor level so the comparison will be Naive vs Sotorasib 
# annotate_data <- column_to_rownames(as.data.frame(annotate_data),'sample') 
# annotate_data <- annotate_data %>% filter(!(condition == "Gemcitabine" & replicate == 2))
# annotation_data$condition <- factor(annotation_data$condition, levels = c("naive", "Gemcitabine", "Sotorasib", "Cisplatin"))
annotation_data$condition <- factor(annotation_data$condition, levels = c("Naive", "Sotorasib"))
write_csv(annotation_data, paste0(output,"annotationmatrix.csv"))

# Start deseq 
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
contrast_info <- c("condition", "Sotorasib", "Naive")
padj_threshold <- 0.05

# Relevel the condition factor to indicate that the reference is naive
dds$condition <- relevel(dds$condition, ref = "Naive")

# Perform DESeq analysis
dds <- DESeq(dds)

# Get results for the specified comparison
full_results <- results(dds, contrast = contrast_info)
saveRDS(full_results, paste0(output, "DESeqResultnvS.rds"))


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
proteinatlas <- read_tsv("/projects/b1042/GoyalLab/aleona/bulk/Protein_data/proteinatlas.tsv", 
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
write_csv(all_results, paste0(output, "DESeqwithgeneid_NaivevsSotfull.csv"))
write_csv(saveabove2, paste0(output, "DESeq_NvSabsLFCresult.csv"))
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

write_csv(over_all, paste0(output, "overexpress.csv"))
write_csv(under_all, paste0(output, "underexpress.csv"))

#Filter just surfaceproteins 
over_sprotein <- over_all %>% filter(`Protein Class`=="Surface Proteins")
under_sprotein <- under_all %>% filter(`Protein Class`=="Surface Proteins")

write_csv(over_sprotein, paste0(output, "overexpress_surfaceprot.csv"))
write_csv(under_sprotein, paste0(output, "underexpress_surfaceprot.csv"))


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
