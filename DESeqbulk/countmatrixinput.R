#Code to extract Count Input Matrix from STAR - created by Aurelia on 2024-08-09 and is based on Keerthana's code 
# last modified on 
library(Seurat)
library(dplyr)
library(tidyverse)


##Main----------------------------------
pathToFiles = "/projects/b1042/GoyalLab/aleona/STARAlignmentforBulk/extData2/NK_H23/results/"
# pathToFiles = "/projects/b1042/GoyalLab/aleona/DESEQ/rawData/bulkRNASeq/"
extData <- "/projects/b1042/GoyalLab/aleona/bulk/extractedData/NK_H23/matrix/"
# extData <- "/projects/b1042/GoyalLab/aleona/DESEQ/extData/rawdata1/"
Scripts<- "/projects/b1042/GoyalLab/aleona/LabPipelines/DESeqbulk /"

source(paste0(Scripts, "deseqbulkfunction.R"))

create_dir_if_not_exists(extData)

#Search all subdirectories uniformly and don't need to distinguish between different levels of subdirectories.
files <- list.files(path = pathToFiles, 
                    pattern = "*ReadsPerGene.out.tab$", 
                    recursive = TRUE, 
                    full.names = TRUE)

# only check for L001
# filtered_files <- files[grepl("L001", files, fixed = TRUE)]
filtered_files <- files

# Define your patterns
# .*?: Matches any character (except newline) zero or more times, as few times as possible (non-greedy).
# (\\w+[-_]\\w+[-_]\\w+): captures group that matches the sample name 
# \\w+: Matches one or more word characters
# [-_]: Matches either a hyphen or an underscore.
# .*?: Matches any remaining characters up to the lane identifier.
match_pattern <- ".*?(\\w+[-_]\\w+[-_]S\\w+).*?"
#Change match pattern as needed 


# Process files
datasets <- map_chr(filtered_files, function(file) {
    extract_sample_name(file, match_pattern)
  }
)

# Extract sample information. Can use Map_chr as well. 
#datasets <- map_chr(basename(files), ~ extract_sample_info(., match_pattern, remove_prefix, remove_suffix))

# Display results
print(datasets)

# [1] "ID_NK_H23-E10-7_cis1_S34"   NA                          
# [3] "ID_NK_H23-E10-7_cis2_S35"   NA                          
# [5] "ID_NK_H23-E10-7_cis3_S40"   NA                          
# [7] "ID_NK_H23-E10-7_Gem1_S33"   NA                          
# [9] "ID_NK_H23-E10-7_Gem2_S39"   NA                          
# [11] "ID_NK_H23-E10-7_Gem3_S41"   NA                          
# [13] "ID_NK_H23-E10-7_naive1_S38" NA                          
# [15] "ID_NK_H23-E10-7_naive2_S42" NA                          
# [17] "ID_NK_H23-E10-7_naive3_S44" NA                          
# [19] "ID_NK_H23-E10-7_sot1_S36"   NA                          
# [21] "ID_NK_H23-E10-7_sot2_S37"   NA                          
# [23] "ID_NK_H23-E10-7_sot3_S43"   NA               

columnnames <- c()
for (i in datasets) {
    columnnames <- c(columnnames,
                     paste0("gene_", i),
                     paste0(i, "_unstranded"),
                     paste0(i, "_forwardstrand"),
                     paste0(i, "_reversestrand")
    )
}

finaltable <- filtered_files %>%
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


write_csv(finaltable, paste0(extData,"finaltable.csv"))
write_csv(strandWiseCounts, paste0(extData,"strandWiseCounts.csv"))
write_csv(countmatrix, paste0(extData,"countmatrix.csv"))


#Adding gene_id to ENSEMBL only data and combine--------------------------------
ensembl.ids <- countmatrix['ensembl_gene']

ensembl_98 <- useEnsembl(biomart = 'genes', 
                         dataset = 'hsapiens_gene_ensembl',
                         version = 98)

gene_id <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                 filters = "ensembl_gene_id",
                 values = ensembl.ids$ensembl_gene,
                 mart = ensembl_98)


all_results <- add_gene_names(countmatrix, gene_id)

write_csv(all_results, paste0(extData,"countmatrixwithgenename.csv"))





          