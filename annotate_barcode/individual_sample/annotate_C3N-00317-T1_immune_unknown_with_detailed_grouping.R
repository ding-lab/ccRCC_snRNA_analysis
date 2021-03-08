# Yige Wu @WashU Apr 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input seurat object
srat1 <- readRDS(file = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/subset_recluster/subset_C3N-00317-T1_immune_unknown_and_recluster/20210305.v1/C3N-00317-T1.Immune_Unknown.Reclustered.20210305.v1.RDS")
srat2 <- readRDS(file = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/subset_recluster/subset_C3N-00317-T1_lymphoid_cells/20210305.v1/C3N-00317-T1.Immune_Unknown.Reclustered.Res2.RDS")

# annotate meta data ----------------------------------------------------------
metadata_df1 <- srat1@meta.data
metadata_df1$barcode <- rownames(metadata_df1)
metadata_df2 <- srat2@meta.data
metadata_df2$barcode <- rownames(metadata_df2)
## combine
metadata_subset_df1 <- metadata_df1 %>%
  filter(!(barcode %in% metadata_df2$barcode))
metadata_subset_df1 <- metadata_df1 %>%
  filter(!(barcode %in% metadata_df2$barcode))