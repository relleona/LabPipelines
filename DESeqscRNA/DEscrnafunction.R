#This is the functions used in all the process of DESEQ2 for scRNA sequencing 
library(Seurat)
library(dplyr)
library(dbplyr)
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(readr)
library(Matrix)
library(patchwork)
library(ggrepel)
library(tidyr)
library(stringr)
library(VennDiagram)
library(gridExtra)

# Load required libraries
required_packages <- c(
  "Seurat", "dbplyr", "dplyr", "Matrix", "tidyverse", "DESeq2", "VennDiagram",
  "biomaRt", "readr", "patchwork", "ggrepel", "tidyr","stringr","gridExtra")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)



# Add this function to perform DESeq2 analysis
perform_seuratDE <- function(seurat_object, condition1, condition2, id_column = NULL) {
  # Make sure that the conditions are in the right metadata variable name  
  if (!is.null(id_column)) {
    # Print available identities
    print("Available identities:")
    print(table(seurat_object[[id_column]]))
  }
  
  # Set the active identity if id_column is provided
  if (!is.null(id_column)) {
    Idents(seurat_object) <- id_column
  }
  
  de.markers <- FindMarkers(seurat_object, ident.1 = condition1, ident.2 = condition2)
  
  return(de.markers)
}



# Add this function to perform DESeq2 analysis -> according to seurat this should be done after pseudobulking 
perform_deseq2_analysis <- function(seurat_object, groups, condition1, condition2, id_column = NULL) {
  # Make sure that the conditions are in the right metadata variable name 
  if (!is.null(id_column)) {
    # Print available identities
    print("Available identities:")
    print(table(seurat_object[[id_column]]))
  }
  
  # Perform pseudobulking
  pseudo_dataset <- AggregateExpression(seurat_object, 
                                        assays = "RNA", 
                                        return.seurat = TRUE, 
                                        group.by = groups)
  
  # Print the groups made
  print("Pseudobulk groups:")
  print(table(Idents(pseudo_dataset)))
  
  # Add the original condition1 and condition2 identities under the name "condition" in the new pseudo_dataset 
  if (!is.null(id_column)) {
    # Create a mapping of cell barcodes to conditions
    condition_map <- setNames(as.character(seurat_object[[id_column]]), colnames(seurat_object))
    
    # Extract the condition for each aggregated group
    pseudo_dataset$condition <- sapply(strsplit(colnames(pseudo_dataset), "_"), function(x) {
      unique_conditions <- unique(condition_map[x])
      if (length(unique_conditions) == 1) {
        return(unique_conditions)
      } else {
        return(paste(unique_conditions, collapse = "_"))
      }
    })
  } else {
    # If no id_column is provided, use the current identity
    condition_map <- setNames(as.character(Idents(seurat_object)), colnames(seurat_object))
    
    pseudo_dataset$condition <- sapply(strsplit(colnames(pseudo_dataset), "_"), function(x) {
      unique_conditions <- unique(condition_map[x])
      if (length(unique_conditions) == 1) {
        return(unique_conditions)
      } else {
        return(paste(unique_conditions, collapse = "_"))
      }
    })
  }
  
  
  # Set the active identity to the new condition column
  Idents(pseudo_dataset) <- "condition"
  
  # Ensure the conditions exist in the pseudobulk dataset
  if (!(condition1 %in% unique(Idents(pseudo_dataset))) || !(condition2 %in% unique(Idents(pseudo_dataset)))) {
    stop("One or both conditions not found in the pseudobulk dataset. Available identities: ", 
         paste(unique(Idents(pseudo_dataset)), collapse = ", "))
  }
  
  # Perform differential expression analysis
  de.markers <- FindMarkers(object = pseudo_dataset, 
                            ident.1 = condition1, 
                            ident.2 = condition2,
                            test.use = "DESeq2")
  
  return(de.markers)
}


match_and_classify_proteins <- function(ordered_data, match_column, proteinatlas) {
  # Prepare proteinatlas data
  proteinatlas_prepared <- proteinatlas %>%
    # Select relevant columns
    dplyr::select(Ensembl, Gene, `Gene synonym`, `Protein class`) %>%
    dplyr::mutate(
      # Convert Gene and Gene synonym to uppercase for case-insensitive matching
      Gene = toupper(Gene),
      `Gene synonym` = toupper(`Gene synonym`),
      # Create boolean columns for each protein class
      `Predicted membrane proteins` = grepl("Predicted membrane proteins", `Protein class`),
      `G-protein coupled receptors` = grepl("G-protein coupled receptors", `Protein class`),
      `Transcription factors` = grepl("Transcription factors", `Protein class`), 
      `CD markers` = grepl("CD markers", `Protein class`),
      # Create a character string of gene names and synonyms
      gene_list = purrr::map2_chr(Gene, `Gene synonym`, ~ paste(unique(c(.x, strsplit(.y, ",\\s*")[[1]])), collapse = "; ") %>% 
                                    stringr::str_trim())
    )
  
  # Perform the matching and classification
  result <- ordered_data %>%
    # Convert the match column to uppercase for case-insensitive matching
    dplyr::mutate(!!match_column := toupper(!!sym(match_column))) %>%
    # Process each row individually
    dplyr::rowwise() %>%
    dplyr::mutate(
      # Find matching rows in proteinatlas_prepared
      matched_info = list(proteinatlas_prepared[stringr::str_detect(proteinatlas_prepared$gene_list, fixed(!!sym(match_column))), ])
    ) %>%
    # Return to processing the entire dataframe at once
    dplyr::ungroup() %>%
    dplyr::mutate(
      # Determine if each protein class is present in the matches
      `Predicted membrane proteins` = purrr::map_lgl(matched_info, ~ any(.$`Predicted membrane proteins`)),
      `G-protein coupled receptors` = purrr::map_lgl(matched_info, ~ any(.$`G-protein coupled receptors`)),
      `Transcription factors` = purrr::map_lgl(matched_info, ~ any(.$`Transcription factors`)),
      `CD markers` = purrr::map_lgl(matched_info, ~ any(.$`CD markers`)),
      # Combine all matched gene names and synonyms
      matched_names = purrr::map_chr(matched_info, ~ paste(.$gene_list, collapse = "; "))
    ) %>%
    dplyr::mutate(
      # Classify proteins based on their characteristics
      `Protein Class` = dplyr::case_when(
        (`Predicted membrane proteins` | `G-protein coupled receptors` | `CD markers`) & `Transcription factors` ~ "Surface Proteins, Transcription Factors",
        `Predicted membrane proteins` | `G-protein coupled receptors` | `CD markers` ~ "Surface Proteins",
        `Transcription factors` ~ "Transcription Factors",
        TRUE ~ "Others"
      )
    ) %>%
    # Remove the temporary matched_info column
    dplyr::select(-matched_info)
  
  # Return the result
  return(result)
}




match_sc_bulk <- function(bulkdata, scrnadata, gene_bulk = "gene_column", gene_scrna = "gene_list") {
  
  # Check if required columns exist in bulkdata
  required_cols <- c(gene_bulk, "Predicted membrane proteins", "G-protein coupled receptors", 
                     "Transcription factors", "CD markers", "Protein Class")
  missing_cols <- setdiff(required_cols, names(bulkdata))
  # Check if gene_scrna column exists in scrnadata
  
  if (!gene_scrna %in% names(scrnadata)) {
    stop(paste("Column", gene_scrna, "not found in scrnadata"))
  }
  
  # Process bulkdata with only the relevant information 
  process_data <- bulkdata[, required_cols]
  
  # Convert gene names to uppercase for case-insensitive matching
  process_data[[gene_bulk]] <- toupper(process_data[[gene_bulk]])
  scrnadata[[gene_scrna]] <- toupper(scrnadata[[gene_scrna]])
  
  result <- process_data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      matched_names = paste(
        scrnadata[[gene_scrna]][stringr::str_detect(scrnadata[[gene_scrna]], fixed(!!sym(gene_bulk)))],
        collapse = "; "
      )
    ) %>%
    dplyr::ungroup()
  
  return(result)
}


summarize_gene_overlap <- function(bulkdata, scrnadata, 
                                   gene_bulk = "external_gene_name", gene_scrna = "matched_names", 
                                   bulk_lfc_cutoff = 2, bulk_p_adj_cutoff = 0.05,
                                   scrna_lfc_cutoff = 2, scrna_p_adj_cutoff = 0.05,
                                   output_folder) {
  
  # Generate dynamic output file names
  output_file <- file.path(output_folder, sprintf("gene_output_summary_lfcbulk%.1f_lfcsc%.1f.txt", bulk_lfc_cutoff, scrna_lfc_cutoff))
  venn_file <- file.path(output_folder, sprintf("gene_overlap_venn_lfcbulk%.1f_lfcsc%.1f.png", bulk_lfc_cutoff, scrna_lfc_cutoff))
  
  # Filter bulkdata based on log2FoldChange and adjusted p-value
  bulkdata_filtered <- bulkdata %>%
    dplyr::filter(abs(log2FoldChange) > bulk_lfc_cutoff) %>%
    dplyr::filter(padj < bulk_p_adj_cutoff)
  
  # Filter scrnadata based on avg_log2FC and p_val_adj
  scrnadata_filtered <- scrnadata %>%
    dplyr::filter(abs(avg_log2FC) > scrna_lfc_cutoff) %>%
    dplyr::filter(p_val_adj < scrna_p_adj_cutoff)
  
  # Use match_sc_bulk to get matched genes
  matched_genes <- match_sc_bulk(bulkdata_filtered, scrnadata_filtered, gene_bulk = gene_bulk, gene_scrna = gene_scrna)
  
  # Extract unique genes from both datasets
  bulk_genes <- unique(toupper(bulkdata_filtered[[gene_bulk]]))
  scrna_genes <- unique(unlist(strsplit(scrnadata_filtered[[gene_scrna]], "; ")))
  scrna_genes <- toupper(scrna_genes[scrna_genes != "NA"])
  
  # Find overlapping genes
  overlapping_genes <- intersect(bulk_genes, scrna_genes)
  
  # Calculate overlap statistics
  total_bulk_genes <- length(bulk_genes)
  total_scrna_genes <- length(scrna_genes)
  overlap_count <- length(overlapping_genes)
  percent_overlap <- (overlap_count / min(total_bulk_genes, total_scrna_genes)) * 100
  
  # Create Venn diagram with outline colors
  venn_plot <- venn.diagram(
    x = list(Bulk = bulk_genes, scRNA = scrna_genes),
    filename = NULL,
    fill = c("lightblue", "lightpink"),
    col = c("blue", "red"),
    alpha = 0.5,
    lwd = 2,
    lty = "solid",
    label.col = "black",
    cex = 1.5,
    fontfamily = "sans",
    cat.col = c("blue", "red"),
    cat.cex = 1.2,
    cat.fontfamily = "sans"
  )
  
  # Save Venn diagram as a PNG file
  png(venn_file, width = 800, height = 800)
  grid.draw(venn_plot)
  dev.off()
  
  # Prepare summary text
  summary_text <- sprintf(
    "Gene Overlap Summary:\n\n" +
      "Total genes in filtered bulk data: %d\n" +
      "Total genes in filtered scRNA data: %d\n" +
      "Number of overlapping genes: %d\n" +
      "Percent overlap: %.2f%%\n" +
      "Bulk LFC cutoff: %.2f\n" +
      "Bulk p-value adjusted cutoff: %.3f\n" +
      "scRNA LFC cutoff: %.2f\n" +
      "scRNA p-value adjusted cutoff: %.3f\n",
    total_bulk_genes, total_scrna_genes, overlap_count, percent_overlap, 
    bulk_lfc_cutoff, bulk_p_adj_cutoff, scrna_lfc_cutoff, scrna_p_adj_cutoff
  )
  
  # Write summary to text file
  writeLines(summary_text, output_file)
  
  # Print summary to console
  cat(summary_text)
  
  # Return a list containing the requested information
  return(list(
    matched_genes = matched_genes,
    summary_file = output_file,
    venn_file = venn_file,
    overlapping_genes = overlapping_genes
  ))
}
