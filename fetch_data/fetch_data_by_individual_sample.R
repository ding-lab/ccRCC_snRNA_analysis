# Yige Wu @WashU Apr 2020
## for plotting the cluster number from the integrated data to individual sample clusters

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
## input path to individual srat objects
srat_paths <- fread(input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20200717.v1/Seurat_Object_Paths.20200717.v1.tsv", data.table = F)

# plot cell types to umap by each aliquot ---------------------------------
metadata_df <- NULL
for (snRNA_aliquot_id_tmp in unique(srat_paths$Aliquot)) {
  ## input seurat object
  srat_path <- srat_paths$Path_box_seurat_object[srat_paths$Aliquot == snRNA_aliquot_id_tmp]
  srat <- readRDS(file = srat_path)
  ## fetch data
  metadata_tmp <- FetchData(srat, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2", "nFeature_RNA", "nCount_RNA", "percent.mito"))
  metadata_tmp$individual_barcode <- rownames(metadata_tmp)
  
  ## merge with super table
  metadata_df <- rbind(metadata_tmp, metadata_df)
}

# write output ------------------------------------------------------------
metadata_df <- metadata_df %>%
  rename(aliquot = orig.ident) %>%
  rename(seurat_cluster_id = ident)
file2write <- paste0(dir_out, "Barcode2MetaData.", run_id, ".tsv")
write.table(x = metadata_df, file = file2write, quote = F, sep = "\t", row.names = F)
