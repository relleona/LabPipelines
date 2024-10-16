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
extract_sample_name <- function(filename, pattern = NULL) {
  # If no pattern is provided, use the default pattern
  if (is.null(pattern)) {
    # .*?: Matches any character (except newline) zero or more times, as few times as possible (non-greedy).
    # (\\w+[-_]\\w+[-_]\\w+): captures group that matches the sample name 
    # \\w+: Matches one or more word characters
    # [-_]: Matches either a hyphen or an underscore.
    # .*?: Matches any remaining characters up to the lane identifier.
    pattern <- ".*?(\\w+[-_]\\w+[-_]\\w+).*?"
  }
  
  # Combine the pattern with the lane identifier
  # full_pattern <- paste0(pattern, lane_identifier)
  
  # Extract the match
  match <- regexpr(pattern, filename, perl = TRUE)
  
  if (match == -1) {
    warning("Pattern not found in the filename.")
    return(NA)
  }
  
  # Extract the captured group
  start <- attr(match, "capture.start") # returns position of the pattern
  length <- attr(match, "capture.length") # length of pattern 
  sample_name <- substr(filename, start, start + length - 1)
  
  return(sample_name)
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



create_annotation_matrix <- function(sample_names, 
                                     sample_type = "paired-end",
                                     condition_mapping = NULL) {
  
  # Create the initial annotation matrix
  annotation_matrix <- data.frame(
    sample = sample_names,
    type = sample_type,
    stringsAsFactors = FALSE
  )
  # Explanation of the regex for replicate extraction:
  # .*          : Matches any characters at the beginning of the string.
  # (?:...)     : This is a non-capturing group. It groups parts of the regex without creating a backreference.
  # [_-]        : Matches either an underscore or a dash.
  # ([A-Za-z])  : Captures any single letter (uppercase or lowercase) in the first capturing group.
  # (\\d)       : Captures a single digit in the second capturing group.
  # "|"         : This is an OR operator. It separates two alternative patterns.
  # [_-](\\d)[_-] : This is the second alternative. It matches an underscore or dash, followed by a digit (captured in the third group), followed by another underscore or dash.
  # .*          : Matches any characters at the end of the string.
  # "\\2\\3"    : In the replacement part, this refers to either the second capture group (the digit after a letter) or the third capture group (the digit between underscores/dashes), whichever was matched.
  # perl = TRUE : This enables Perl-compatible regular expressions, which support more advanced features.
  
  # Function to extract replicate from sample name
  extract_replicate <- function(sample) {
    as.integer(sub(".*(?:[_-]([A-Za-z])(\\d)[_-]|[_-](\\d)[_-]).*", "\\2\\3", sample, perl = TRUE))
  }
  
  # Function to extract condition from sample name
  extract_condition <- function(sample) {
    sub("^(.+?)_.*", "\\1", sample)
  }
  
  if (is.null(condition_mapping)) {
    # Case 1: No condition mapping provided
    annotation_matrix <- annotation_matrix %>%
      mutate(
        condition = sapply(sample, extract_condition),
        replicate = sapply(sample, extract_replicate)
      )
  } else if (is.character(condition_mapping)) {
    # Case 2: Condition mapping with only conditions
    annotation_matrix <- annotation_matrix %>%
      mutate(
        condition = condition_mapping[sample],
        replicate = sapply(sample, extract_replicate)
      )
  } else if (is.list(condition_mapping)) {
    # Case 3: Condition mapping with both conditions and replicates
    annotation_matrix <- annotation_matrix %>%
      mutate(
        condition = sapply(sample, function(s) condition_mapping[[s]]$condition),
        replicate = sapply(sample, function(s) condition_mapping[[s]]$replicate)
      )
  } else {
    stop("Invalid condition_mapping format")
  }
  
  # Set row names
  rownames(annotation_matrix) <- annotation_matrix$sample
  
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


# Just adding the gene names 
add_gene_names <- function(count_matrix, gene_id_df) {
  # Load required package
  # if (!requireNamespace("dplyr", quietly = TRUE)) {
  #   stop("dplyr package is required but not installed. Please install it using install.packages('dplyr')")
  # }
  # 
  # Ensure count_matrix has ensembl_gene as a column
  if (!("ensembl_gene" %in% colnames(count_matrix))) {
    count_matrix <- dplyr::mutate(count_matrix, ensembl_gene = rownames(count_matrix))
  }
  
  # Trim whitespace and convert to uppercase for both dataframes
  count_matrix$ensembl_gene <- toupper(trimws(count_matrix$ensembl_gene))
  gene_id_df$ensembl_gene_id <- toupper(trimws(gene_id_df$ensembl_gene_id))
  
  # Perform the join
  result <- count_matrix %>%
    dplyr::left_join(gene_id_df, by = c("ensembl_gene" = "ensembl_gene_id")) %>%
    dplyr::select(ensembl_gene, external_gene_name, dplyr::everything())
  
  # Count and print the number of NAs
  na_count <- sum(is.na(result$external_gene_name))
  total_count <- nrow(result)
  print(paste("Number of NAs:", na_count, "out of", total_count, "total rows"))
  
  # Print some examples of unmatched gene IDs
  if (na_count > 0) {
    unmatched <- result %>% 
      dplyr::filter(is.na(external_gene_name)) %>% 
      dplyr::select(ensembl_gene) %>% 
      head(5)
    print("Examples of unmatched gene IDs:")
    print(unmatched)
  }
  
  return(result)
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


