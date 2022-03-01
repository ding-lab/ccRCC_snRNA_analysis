# Yige Wu @WashU Nov 2020

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
## input barcode2UMAP
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/recluster/recluster_tumorcells/downsample_fixednumber_and_recluster_tumor_cells_in_selected_samples_katmai/20220222.v1/UMAPData.2000TumorCellReclustered.20220222.v1.tsv")
## input meta data
### accidently use the wrong sample id
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# count -------------------------------------------------------------------
barcode2umap_df$easy_id <- mapvalues(x = barcode2umap_df$orig.ident, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Aliquot.snRNA.WU))
cellnumber_percluster_df <- barcode2umap_df %>%
  mutate(id_cluster_uniq = paste0(easy_id, "_C", seurat_clusters)) %>%
  select(id_cluster_uniq) %>%
  table() %>%
  as.data.frame() %>%
  rename(id_cluster_uniq = '.') %>%
  mutate(easy_id = str_split_fixed(string = id_cluster_uniq, pattern = "_C", n = 2)[,1]) %>%
  arrange(easy_id)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellNumberPerTumorSeuratCluster.", run_id, ".tsv")
write.table(x = cellnumber_percluster_df, file = file2write, quote = F, sep = "\t", row.names = F)
