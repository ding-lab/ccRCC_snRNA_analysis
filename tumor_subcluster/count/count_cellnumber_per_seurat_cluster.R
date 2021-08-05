# Yige Wu @WashU Nov 2020

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
## input barcode2UMAP
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_individual_tumorcellreclustered_on_katmai/20201127.v1/MetaData_TumorCellOnlyReclustered.20201127.v1.tsv")

# count -------------------------------------------------------------------
cellnumber_percluster_df <- barcode2umap_df %>%
  mutate(id_cluster_uniq = paste0(easy_id, "_C", seurat_clusters)) %>%
  select(id_cluster_uniq) %>%
  table() %>%
  as.data.frame() %>%
  rename(id_cluster_uniq = '.') %>%
  mutate(easy_id = str_split_fixed(string = id_cluster_uniq, pattern = "_C", n = 2)[,1])

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellNumberPerTumorSeuratCluster.", run_id, ".tsv")
write.table(x = cellnumber_percluster_df, file = file2write, quote = F, sep = "\t", row.names = F)
