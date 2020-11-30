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
## input seurat object paths
srat_paths <- fread(input = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)

# process by each aliquot ----------------------------------------------------
barcode_metadata_df <- NULL
for (aliquot_tmp in srat_paths$Aliquot) {
  ## input srat object
  srat_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  srat <- readRDS(file = srat_path)
  
  ## extract current meta data
  barcode_metadata_tmp <- FetchData(object = srat, vars = c("orig.ident", "seurat_clusters", "UMAP_1", "UMAP_2"))
  barcode_metadata_tmp$barcode <- rownames(barcode_metadata_tmp)
  ## bind with the super table
  barcode_metadata_df <- rbind(barcode_metadata_tmp, barcode_metadata_df)
}
## reformat
barcode_metadata_df <- barcode_metadata_df %>%
  rename(id_seurat_cluster = seurat_clusters)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Barcode2SeuratClusterID.", run_id, ".tsv")
write.table(x = barcode_metadata_df, file = file2write, sep = '\t', quote = F, row.names = F)
  