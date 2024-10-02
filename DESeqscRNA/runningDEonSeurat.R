#Define paths
main_path <- "/projects/b1042/GoyalLab/aleona/DESeqscRNA/"
extData <- paste0(main_path, "extractedData/")
rawData <- paste0(main_path, "rawData/")
scripts <- paste0(main_path, "Scripts/")
plotpath <- paste0(main_path, "plot/")

# #Run the source function scripts
source(paste0(scripts, "DEscrnafunction.R"))
source(paste0(scripts, "seuratpipelinefunctions.R"))

#Create directory
create_dir_if_not_exists(extData)
create_dir_if_not_exists(rawData)
create_dir_if_not_exists(scripts)
create_dir_if_not_exists(plotpath)

#Input the integrated dataset 
integrated_seurat <- readRDS(paste0(extData, "integrated_naiveSot.rds"))

##Differential Gene Expression Analysis--------------------------
#Using normal FindMarkers 
dge_results <- perform_seuratDE(integrated_seurat, "Sot_078", "Naive_078", id_column = "orig.ident")

print(head(dge_results))

# #DESEQ  -> pseudobulk 
# deseq2_results <- perform_deseq2_analysis(integrated_seurat, groups= c("orig.ident", "seurat_clusters"), "Sot_078", "Naive_078", id_column = "orig.ident")
# 
# # Display top differentially expressed genes
# print(head(deseq2_results$de.markers))


##Match Protein Classes ------------------------------------------

#Find the ensembl gene that corresponds to the rows of the dge_results 

#Read the Protein atlast database 
proteinatlas <- read_tsv("/projects/b1042/GoyalLab/aleona/DESeqbulk/Protein_data/proteinatlas.tsv", 
                         col_names = TRUE, 
                         show_col_types = FALSE)

#Take out the rownames and add it to a column called Gene 
results <- rownames_to_column(dge_results, "gene_column")

# Assuming proteinatlas is already loaded
full_results <- match_and_classify_proteins(results, "gene_column", proteinatlas)

# Convert to another dataset
ordered <- as.data.frame(full_results)

# Convert to absolute values
ordered$avg_log2FC <- abs(ordered$avg_log2FC)

# Order the data frame
ordered <- ordered[order(ordered$avg_log2FC, decreasing = TRUE), ]

#Only the Log2FoldChange above 2 
saveabove2 <- ordered %>% filter(avg_log2FC > 2 ) %>% filter(p_val_adj < 0.05 )

#Saving results 
# write_csv(saveabove2, paste0(extData2, "keeptrackSep32024.csv"))
write_csv(results, paste0(extData, "DESeqwithgeneid_NaivevsSotfullscrna.csv"))
write_csv(saveabove2, paste0(extData, "DESeq_NvSabsLFCresultscrna.csv"))

#Filter the data to over and under expressed------------------------------------------
ordered1 <- as.data.frame(full_results)

#Split between overexpression and underexpression 
# Over is positive because we want that is more highly expressed in Sotorasib 
over <- ordered1 %>% filter(avg_log2FC > 2 ) %>% filter(p_val_adj < 0.05 ) 
under <- ordered1 %>% filter(avg_log2FC < -2 ) %>% filter(p_val_adj < 0.05 )

write_csv(over, paste0(extData, "overexpress.csv"))
write_csv(over, paste0(extData, "underexpress.csv"))

#Filter just surfaceproteins 
over_sprotein <- over %>% filter(`Protein Class`=="Surface Proteins")
under_sprotein <- under %>% filter(`Protein Class`=="Surface Proteins")

write_csv(over_sprotein, paste0(extData, "overexpress_surfaceprot.csv"))
write_csv(under_sprotein, paste0(extData, "underexpress_surfaceprot.csv"))


#Compare the dataset with bulk and find overlaps---------------------
#Bulk dataset initialize the path 
bulk_path <- "/projects/b1042/GoyalLab/aleona/DESeqbulk/extData/rawData2/"

##Match it with bulk data-overexpress-----------------
bulkover_data <- read_csv(paste0(bulk_path,"overexpress_surfaceprot.csv"))
bulkover_data <- as.data.frame(bulkover_data)

matched_genes <- match_sc_bulk(bulkover_data, over_sprotein, gene_bulk = "external_gene_name", gene_scrna = "matched_names")

write_csv(matched_genes, paste0(extData, "match_sc_bulk_overexpress.csv"))

##Match it with bulk data-underexpress 
bulkunder_data <- read_csv(paste0(bulk_path, "underexpress_surfaceprot.csv"))
bulkunder_data <- as.data.frame(bulkunder_data)

matched_genes <- match_sc_bulk(bulkunder_data, under_sprotein, gene_bulk = "external_gene_name", gene_scrna = "matched_names")

write_csv(matched_genes, paste0(extData, "match_sc_bulk_underexpress.csv"))


##Match it with bulk data
gene_overlaps <- paste0(extData, "bulkvscrna/")
create_dir_if_not_exists(gene_overlaps)

bulk_data <- read_csv(paste0(bulk_path, "DESeqwithgeneid_NaivevsSotfull2.csv"))
bulk_data <- as.data.frame(bulk_data)

scrna_data <- as.data.frame(full_results)

result1 <- summarize_gene_overlap(bulkdata= bulk_data, scrnadata=scrna_data , gene_bulk = "external_gene_name", gene_scrna = "matched_names", 
                       bulk_lfc_cutoff = 2, bulk_p_adj_cutoff = 0.05,
                       scrna_lfc_cutoff = 2, scrna_p_adj_cutoff = 0.05,
                       output_folder=gene_overlaps)



















