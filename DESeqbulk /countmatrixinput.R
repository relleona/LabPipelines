#Code to extract Count Input Matrix from STAR - created by Aurelia on 2024-08-09 and is based on Keerthana's code 
# last modified on 
library(Seurat)
library(dplyr)
library(tidyverse)


##Main----------------------------------
pathToFiles = "/projects/b1042/GoyalLab/aleona/STARAlignmentPipeline/extData4/bulkRNAseq/results/"
# pathToFiles = "/projects/b1042/GoyalLab/aleona/DESEQ/rawData/bulkRNASeq/"
extData <- "/projects/b1042/GoyalLab/aleona/DESEQ/extData/bulkRNASeq2/"
# extData <- "/projects/b1042/GoyalLab/aleona/DESEQ/extData/rawdata1/"
Scripts<- "/projects/b1042/GoyalLab/aleona/DESEQ/Script/"

source(paste0(Scripts, "deseqbulkfunction.R"))

#Search all subdirectories uniformly and don't need to distinguish between different levels of subdirectories.
files <- list.files(path = pathToFiles, 
                    pattern = "*ReadsPerGene.out.tab$", 
                    recursive = TRUE, 
                    full.names = TRUE)

# Define your patterns
match_pattern <- "_(.+?)_(\\d+)_(S\\d+)"
remove_prefix <- "(^.*_E7)"  # Anything before E7
remove_suffix <- "_L001ReadsPerGene\\.out\\.tab"
remove_prefix2 <- "(^.*_MP2)"  # Anything before MP2

# Process files
datasets <- map_chr(files, function(file) {
  if (str_detect(file, remove_prefix)) {
    extract_sample_info(file, match_pattern, remove_prefix, remove_suffix)
  } else if (str_detect(file, remove_prefix2)) {
    extract_sample_info(file, match_pattern, remove_prefix2, remove_suffix)
  } else {
    extract_sample_info(file, match_pattern)
  }
})

# Extract sample information. Can use Map_chr as well. 
#datasets <- map_chr(basename(files), ~ extract_sample_info(., match_pattern, remove_prefix, remove_suffix))

# Display results
print(datasets)

columnnames <- c()
for (i in datasets) {
  columnnames <- c(columnnames,
                   paste0("gene_", i),
                   paste0(i, "_unstranded"),
                   paste0(i, "_forwardstrand"),
                   paste0(i, "_reversestrand")
  )
}


finaltable <- files %>%
  # applies read_tsv function for each path in files 
  # Skip the first 4 lines and not use the first non_skipped lines as column names 
  map(read_tsv,  skip = 4, col_names = FALSE ) %>%
  # reduce() combines a list of objects into one 
  # cbind is the function used by reduce to combine column_wise, this is under the assumption rows in each file correspond to the same genes and are in the same order. 
  # This is typically the case for RNA-seq count data generated from the same reference genome and annotation.
  purrr::reduce(cbind) 

names(finaltable) <- columnnames
rm(columnnames, datasets)

finaltable <- finaltable %>%
  # Ensure the dataframe has column names
  setNames(c("ensembl_gene", names(.)[-1])) %>%
  # Remove columns starting with "gene"
  dplyr::select(-starts_with("gene"))

finaltable %>% head()

#'The purpose of this calculation is to assess the strand specificity of your RNA-seq library preparation method:
#'- In unstranded RNA-seq, reads align to both strands approximately equally. 
#'- In stranded RNA-seq, reads should align predominantly to either the forward or reverse strand, depending on the library preparation method.
#'
#'The propF and propR calculations give you the proportion of reads aligning to the forward and reverse strands, respectively, compared to the unstranded count.
#'Interpreting the results:
#'- For unstranded data, you'd expect propF and propR to be close to 0.5 each.
#'- For forward-stranded data, propF should be close to 1 and propR close to 0.
#'- For reverse-stranded data, propF should be close to 0 and propR close to 1.
#'
strandWiseCounts <- finaltable %>%
  # Remove the ensembl_gene column
  dplyr::select(-ensembl_gene) %>%  # Use dplyr::select to ensure the right function is used
  # Sum up counts for each sample across all genes
  summarise(across(everything(), sum)) %>%
  # Reshape data to long format: one row per sample
  pivot_longer(cols = everything(), names_to = "library", values_to = "counts") %>%
  # Separate the library name into two parts: 
  # 'combined' (everything up to the last underscore) and 'stranding'
  separate(col = library, into = c("combined", "stranding"), sep = "_(?=[^_]+$)", extra = "merge") %>%
  # Further separate 'combined' into 'treatment' and 'replicate'
  separate(col = combined, into = c("treatment", "replicate"), sep = "_", extra = "merge") %>%
  # Reshape data to have one column for each strand type (unstranded, forwardstrand, reversestrand)
  pivot_wider(names_from = stranding, values_from = counts) %>%
  # Calculate proportions for forward and reverse strands
  mutate(
    propF = forwardstrand / unstranded,  # Calculate the proportion of forwardstrand to unstranded
    propR = reversestrand / unstranded   # Calculate the proportion of reversestrand to unstranded
  )

#Create Annotation matrix
countmatrix <- finaltable

#It is found from the StrandWiseCounts data that we should use the unstranded data 
countmatrix <- countmatrix %>% 
  dplyr::select(-ends_with("_reversestrand"), -ends_with("forwardstrand")) 

#Remove the "_unstranded" name 
names(countmatrix) <- str_replace(names(countmatrix), "_unstranded", "")


create_dir_if_not_exists(extData)
write_csv(finaltable, paste0(extData,"finaltable.csv"))
write_csv(strandWiseCounts, paste0(extData,"strandWiseCounts.csv"))
write_csv(countmatrix, paste0(extData,"countmatrix.csv"))


          