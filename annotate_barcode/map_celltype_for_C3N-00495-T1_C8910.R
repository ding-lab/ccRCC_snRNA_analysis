# Yige Wu @WashU Apr 2020
## make barcode to cell type mapping table for the integrated dataset
## just for normal epithelial cells

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(dplyr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------------
## input cluster2celltype
cluster2celltype_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Individual_AllClusters/C3N-00495-T1_C8910_Reclustering.xlsx")
## inpus seurat object
srat <- readRDS(file = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/subset_C3N-00495_C8910_and_recluster/20200922.v1/TumorLikeCells.Reclustered.20200922.v1.RDS")

# map cell type to barcode ------------------------------------------------
barcode2cluster_df <- srat@meta.data 
nrow(barcode2cluster_df)
barcode2cluster_df <- as.data.frame(barcode2cluster_df)
barcode2cluster_df$individual_barcode <- rownames(barcode2cluster_df)
unique(barcode2cluster_df$seurat_clusters)
barcode2cluster_df$seurat_clusters <- as.numeric(as.vector(barcode2cluster_df$seurat_clusters))
unique(barcode2cluster_df$seurat_clusters)
barcode2celltype_df <- merge(barcode2cluster_df, 
                             cluster2celltype_df, 
                             by.x = c("seurat_clusters"), by.y = c("Cluster"), all.x = T)
unique(barcode2cluster_df$seurat_clusters)

## format

# write output ------------------------------------------------------------
nrow(barcode2celltype_df) # 405
file2write <- paste0(dir_out, "Barcode2CellType.", "C3N-00495-T1_C8910.", run_id, ".tsv")
write.table(x = barcode2celltype_df, file = file2write, quote = F, sep = "\t", row.names = F)
  