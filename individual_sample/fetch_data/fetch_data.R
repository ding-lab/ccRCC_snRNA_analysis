# Yige Wu @WashU Apr 2020
## for plotting the cluster number from the integrated data to individual sample clusters

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input path to individual srat objects
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_sample/write_path_to_seurat_objects_on_box/20200219.v1/Seurat_Object_Paths.20200219.v1.tsv", data.table = F)

# plot cell types to umap by each aliquot ---------------------------------
metadata_df <- NULL
for (snRNA_aliquot_id_tmp in unique(srat_paths$Aliquot)) {
  ## input seurat object
  srat_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == snRNA_aliquot_id_tmp]
  srat <- readRDS(file = srat_path)
  
  ## fetch data
  metadata_tmp <- FetchData(srat, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  metadata_tmp$individual_barcode <- rownames(metadata_tmp)
  
  ## merge with super table
  metadata_df <- rbind(metadata_tmp, metadata_df)
}

# write output ------------------------------------------------------------
metadata_df <- metadata_df %>%
  rename(aliquot = orig.ident) %>%
  rename(seurat_cluster_id = ident)
file2write <- paste0(dir_out, "individual_cluster_meta_data.", run_id, ".tsv")
write.table(x = metadata_df, file = file2write, quote = F, sep = "\t", row.names = F)
