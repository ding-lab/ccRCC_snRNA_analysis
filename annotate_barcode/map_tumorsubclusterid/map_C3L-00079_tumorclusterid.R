# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the directory hosting the reclustered object
dir_srats <- "./Resources/snRNA_Processed_Data/Tumor_Subclusters/Reclustered_Objects/"
## input reclustered seurat cluster id to manual cluster id
sratcluster2manualcluster_df <- readxl::read_xlsx(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Tumor_Subcluster/Individual.TumorCluster2Cell_Type.20201008.v1.xlsx")

# process -----------------------------------------------------------------
## get all paths to the seurat objects
paths_srat <- list.files(path = dir_srats, full.names = T)
paths_srat
## process by path
barcodeinfo_all_df <- NULL
for (path_srat in paths_srat) {
  ## input the srat object
  srat <- readRDS(file = path_srat)
  ## extract barcode info
  barcodeinfo_df <- FetchData(object = srat, vars = c("orig.ident", "seurat_clusters", "UMAP_1", "UMAP_2"))
  barcodeinfo_df$barcode_tumorcellreclustered <- rownames(barcodeinfo_df)
  ## comine
  barcodeinfo_all_df <- rbind(barcodeinfo_df, barcodeinfo_all_df)
}

# merge the manual cluster id ---------------------------------------------
# barcodeinfo_all_df$seurat_clusters
# sratcluster2manualcluster_df$id_seurat_cluster
# table(barcodeinfo_all_df$seurat_clusters)
## make the seurat id from factor to vector for merging
barcodeinfo_all_df$seurat_clusters <- as.numeric(as.vector(barcodeinfo_all_df$seurat_clusters))
barcodeinfo_merged_df <- merge(x = barcodeinfo_all_df, y = sratcluster2manualcluster_df,
                            by.x = c("orig.ident", "seurat_clusters"), by.y = c("Aliquot", "id_seurat_cluster"), all.x = T)
## check
table(barcodeinfo_merged_df[, c("seurat_clusters", "id_manual_cluster")])

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "TumorCellReclustered.BarcodeInfo.", run_id, ".tsv")
write.table(x = barcodeinfo_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)


