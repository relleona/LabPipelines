#Script to extract DESeq2 data from scRNAseq data 
#Some functions are taken from Keerthana's code - some are slightly modified

# Load necessary libraries
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis) 
library(RColorBrewer)
library(reticulate)
library(readr)
library(ggrepel)
library(patchwork)
library(DESeq2)
library(Matrix)


#Main--------------------------
#Set the seed so the results are always the same or at least similar 
set.seed(10403)

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

#Load the input data 
naive_raw <- paste0(rawData, "SA2_078_naive/filtered_feature_bc_matrix/")
sot_raw <- paste0(rawData, "SA2_078_SotA/filtered_feature_bc_matrix/")

#Initialize Seuratobject 
#Inputs will be the path to the data and the label you wanted it to be in the Seurat object 
naive <- initSeuratObject(naive_raw, "Naive_078")
SotA <- initSeuratObject(sot_raw, "Sot_078")

##Preprocessing workflow---------------------------------------------
#Make a folder called QCplots ti save the plots created here 
qcpath <- paste0(plotpath, "QCplots/")
create_dir_if_not_exists(qcpath)

#List the paramters that we need to filter the raw data 
# Visualization of the cutoffs will be done in the line below 
# Percentmitocutoff is the maximum mitochondria percentage 
# minCounts and maxCounts are the minimum and maximum 
params <- list(
  naive= list(seuratObject = naive, prefix = "naive", percentMitoCutoff = 10, minCounts = 7000, minFeatures = 2000, maxCounts = Inf, maxFeatures = Inf),
  SotA = list(seuratObject = SotA, prefix = "SotA", percentMitoCutoff = 8, minCounts = 7000, minFeatures = 1500, maxCounts = Inf, maxFeatures = Inf)
)

#Plot all the parameters and look at the min and max range to ensure that we filter the noise 
#This will show if we need to change the parameters above 
plotAll(params, qcpath)

# Apply filters
naive_filtered <- applyFilter(naive, params$naive)
SotA_filtered <- applyFilter(SotA, params$SotA)

#Run the integration workflow 
#Create a folder called integration_plots to save the plots obtained here 
intpath <- paste0(plotpath, "integration_plots/")
create_dir_if_not_exists(intpath)

#Make an object list as an input to the function and run the function 
object_list <- list(naive = naive_filtered, sot = SotA_filtered)
integrated_seurat <- integrate_and_normalize(object_list, n_pcs = 50, intpath)

#Run UMAP functions, 
# n_pcs= indicating the number of Principle components to use
# reduction_name = name of the reduction process that we want to use to plot 
# path_to_folder = indicates which folder to save the umap plot 
final_seurat <- runumap(integrated_seurat, n_pcs = 50, path_to_folder= intpath, reduction_name="harmony", custom_name= "integrated")

#Save the integrated data 
saveRDS(final_seurat, file = paste0(extData, "integrated_naiveSot.rds"))


## Plotting Feature Plot -----------------------------

#Create a folder called the feature_plots to save all the plots created 
fpath <- paste0(plotpath, "feature_plots/")
create_dir_if_not_exists(fpath)

#Plot the entire integrated_seurat
# Dimensional reduction plot, with cells colored by a quantitative feature Defaults to UMAP if available
features = c("PSG9", "ADAM23", "PLAAT4", "APOL1", "ATP10D", "NEO1", "C14orf132", "ELOVL7", "APOL6")

#This will output the feature plot
# object is the seurat object 
# list_fatures is the feature list you wish to plot 
# Path_to_folder is the path to the folder you want the plots to be saved in 
plot_integrated <- runfeatureplot(object= final_seurat, list_features= features , path_to_folder= fpath, custom_name= "integrated")

#if you want to view the plot
print(plot_integrated)

#Split the Seurat object based on the two types 
#Subset or split by the metadata column called orig.ident where it contains labels of both the two different data 
# orig.ident == "Naive_078" means orig.ident is equals to "Naive_078"
#Custom name is what you want to name the plot
naive_int <- subset(final_seurat, subset= orig.ident == "Naive_078")
plot_naive <- runfeatureplot(object= naive_int, list_features= features , path_to_folder= fpath, custom_name= "Naive")

print(plot_naive)

#Do the same for sot
Sot_int<- subset(final_seurat,  subset= orig.ident ==  "Sot_078")
plot_Sot <- runfeatureplot(object= Sot_int, list_features= features , path_to_folder= fpath, custom_name= "Sot")

print(plot_Sot)






