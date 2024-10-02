#This is the functions used in all the process of DESEQ2 for bulkRNA sequencing 
library(Seurat)
library(dplyr)
library(dbplyr)
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(readr)

# Load required libraries
required_packages <- c("Seurat", "dbplyr", "dplyr", "Matrix", "tidyverse", "DESeq2", "biomaRt", "readr")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)


### Functions-----------
#' Extract Sample Information from Filename
#'
#' This function processes a filename to extract relevant sample information
#' based on a given pattern, optionally removing specified prefixes and suffixes.
#'
#' @param filename A character string. The full filename to process.
#'
#' @param match_pattern A regular expression string. The pattern to match and 
#'        extract parts from the cleaned filename. Use capturing groups () to 
#'        specify the parts to extract.
#'
#' @param remove_prefix A character string or NULL. Optional. A regular expression 
#'        pattern to remove from the start of the filename. Default is NULL (no removal).
#'
#' @param remove_suffix A character string or NULL. Optional. A regular expression 
#'        pattern to remove from the end of the filename. Default is NULL (no removal).
#'
#' @return A character string containing the extracted and combined parts of the 
#'         filename based on the match_pattern. Returns NA if the pattern doesn't match.
#'
#' @examples
#' # Basic usage
#' extract_sample_info("SA_MP2_E7_naive_1_S1.tab", "(SA.+?)_(.+?)_(\\d+)_(S\\d+)")
#' # Returns: "SA_MP2_E7_naive_1_S1"
#'
#' # With prefix and suffix removal
#' extract_sample_info("ID:SA_MP2_E7_naive_1_S1_suffix.tab", 
#'                     "(SA.+?)_(.+?)_(\\d+)_(S\\d+)", 
#'                     remove_prefix = "ID:", 
#'                     remove_suffix = "_suffix\\.tab")
#' # Returns: "SA_MP2_E7_naive_1_S1"
#'
#' @note The function will issue a warning if the match_pattern doesn't match 
#'       the cleaned filename.
extract_sample_info <- function(filename, match_pattern, remove_prefix = NULL, remove_suffix = NULL) {
  cleaned_name <- filename
  
  if (!is.null(remove_prefix)) {
    cleaned_name <- str_remove(cleaned_name, paste0("^", remove_prefix))
  }
  
  if (!is.null(remove_suffix)) {
    cleaned_name <- str_remove(cleaned_name, paste0(remove_suffix, "$"))
  }
  
  str_match(cleaned_name, match_pattern)[-1] %>% 
    .[!is.na(.)] %>% 
    paste(collapse = "_")
}


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



###
create_annotation_matrix <- function(sample_names, 
                                     sample_type = "paired-end",
                                     condition_mapping = NULL, patterns) {
  
  # Create the initial annotation matrix
  annotation_matrix <- data.frame(
    sample = sample_names,
    type = sample_type,
    stringsAsFactors = FALSE
  )
  
  # Add condition based on sample names
  if (is.null(condition_mapping)) {
    # Use the first value before underscore as condition
    annotation_matrix <- annotation_matrix %>%
      mutate(condition = sub("^(.+?)_.*", "\\1", sample))
  } else {
    # Use provided condition mapping
    annotation_matrix <- annotation_matrix %>%
      mutate(condition = condition_mapping[sample])
  }
  # Add replicate number
  annotation_matrix <- annotation_matrix %>%
    mutate(replicate = as.integer(sub(".*_(\\d+)_.*", "\\1", sample)))
  
  # Set row names and remove sample column
  rownames(annotation_matrix) <- annotation_matrix$sample
  # annotation_matrix$sample <- NULL
  
  return(annotation_matrix)
}


perform_deseq_analysis <- function(dds, contrast_info, padj_threshold = 0.05) {
  # # Create a copy of the DESeqDataSet
  # dds_copy <- dds
  
  # Relevel the condition factor
  dds$condition <- relevel(dds_copy$condition, ref = reference_condition)
  
  # Perform DESeq analysis
  dds <- DESeq(dds)
  
  # Get results for the specified comparison
  res <- results(dds_copy, contrast = contrast_info)
  
  # Do results function based on adjusted p-value
  res_padj <- results(dds_copy, contrast = contrast_info, alpha = padj_threshold)
  
  # Return both the full results and the padj results
  return(list(full_results = res, padj_results = res_padj))
}


combine_deseq_results <- function(result_list, gene_id_df, comparisons, padj_threshold = 0.05) {
  # Check if the lengths of result_list and comparisons match
  if (length(result_list) != length(comparisons)) {
    stop("The number of results and comparisons must be the same.")
  }
  
  # Process each result
  processed_results <- mapply(function(result, comparison) {
    result %>%
      # Subset DESeqResults object
      # subset(padj < padj_threshold) %>%
      # Convert to data frame and add rownames as a column in one step
      data.frame(ensembl_gene = rownames(.), .) %>%
      # Join with gene_id_df to add more gene information
      # The "ensembl_gene" column is matched with "ensembl_gene_id" in gene_id_df
      left_join(gene_id_df, by = c("ensembl_gene" = "ensembl_gene_id")) %>%
      
      # Add a new column with the comparison name
      mutate(comparison = comparison)
  }, result_list, comparisons, SIMPLIFY = FALSE)
  
  # Combine all processed results into a single dataframe
  # do.call(rbind, ...) is used to efficiently row-bind multiple dataframes
  all_results <- do.call(rbind, processed_results)
  
  return(all_results)
}


match_and_classify_proteins <- function(ordered_data, match_column, proteinatlas) {
  proteinatlas_prepared <- proteinatlas %>%
    dplyr::select(Ensembl, Gene, `Gene synonym`, `Protein class`) %>%
    mutate(
      `Predicted membrane proteins` = grepl("Predicted membrane proteins", `Protein class`),
      `G-protein coupled receptors` = grepl("G-protein coupled receptors", `Protein class`),
      `Transcription factors` = grepl("Transcription factors", `Protein class`), 
      `CD markers` = grepl("CD markers", `Protein class`)
    )
  
  ordered_data %>%
    left_join(proteinatlas_prepared, by = setNames("Ensembl", match_column)) %>%
    mutate(
      `Predicted membrane proteins` = coalesce(`Predicted membrane proteins`, FALSE),
      `G-protein coupled receptors` = coalesce(`G-protein coupled receptors`, FALSE),
      `Transcription factors` = coalesce(`Transcription factors`, FALSE),
      `CD markers` = coalesce(`CD markers`, FALSE),  # Add this line
      `Protein Class` = case_when(
        (`Predicted membrane proteins` | `G-protein coupled receptors` | `CD markers`) & `Transcription factors` ~ "Surface Proteins, Transcription Factors",
        `Predicted membrane proteins` | `G-protein coupled receptors` | `CD markers` ~ "Surface Proteins",
        `Transcription factors` ~ "Transcription Factors",
        TRUE ~ "Others"
      )
    ) %>%
    dplyr::select(-Gene, -`Gene synonym`, -`Protein class`)  # Remove unnecessary columns
}


