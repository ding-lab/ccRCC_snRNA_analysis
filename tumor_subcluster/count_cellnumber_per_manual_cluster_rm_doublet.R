# Yige Wu @WashU Apr 2021

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
## input the barcode-to-tumorsubcluster table
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20201130.v1/Barcode2TumorSubclusterId.20201130.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20200902.v1/scrublet.run20200902_adj_cutoff.united_outputs.tsv", data.table = F)

# count -------------------------------------------------------------------
barcode_merged_df <- merge(x = barcode2tumorsubcluster_df, y = barcode2scrublet_df, by.x = c("orig.ident", "barcode"), by.y = c("Aliquot", "Barcodes"), all.x = T)
cellnumber_percluster_df <- barcode_merged_df %>%
  filter(!predicted_doublet) %>%
  mutate(id_cluster_uniq = Cluster_Name) %>%
  select(id_cluster_uniq) %>%
  table() %>%
  as.data.frame() %>%
  rename(id_cluster_uniq = '.') %>%
  mutate(easy_id = str_split_fixed(string = id_cluster_uniq, pattern = "_C", n = 2)[,1])

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellNumberPerTumorManualCluster.", run_id, ".tsv")
write.table(x = cellnumber_percluster_df, file = file2write, quote = F, sep = "\t", row.names = F)
