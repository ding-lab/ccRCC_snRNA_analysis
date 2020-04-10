# Yige Wu @WashU March 2020
## for each individual sample tumor cell reclustered, plot UMAP

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input gene to cell type table
gene2celltype_df <- fread(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/Gene2CellType_Tab.20200406.v1.tsv", data.table = F)
## input srat object
srat_path <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_integrated_data/clustering/recluster_normal_epithelial_cells_on_local/20200406.v1/normal_epithelial_cells.reclustered.20200406.v1.RDS"
srat <- readRDS(file = srat_path)

# fetch data  ----------------------------------------------------
metadata_df <- srat@meta.data
metadata_df$normal_integrated_barcode <- rownames(metadata_df)

## get the genes within the cell type marker table
## save plot
file2write <- paste0(dir_out, "normal_epithelial_reclustered.", ".metadata.", run_id, ".tsv")
write.table(x = metadata_df, file = file2write, quote = F, sep = "\t", row.names = F)
