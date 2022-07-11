# Yige Wu @WashU Apr 2021

# Yige Wu @WashU Mar 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the barcode-to-tumorsubcluster table
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)

# count -------------------------------------------------------------------
barcode_merged_df <- merge(x = barcode2tumorsubcluster_df, y = barcode2scrublet_df, by.x = c("orig.ident", "barcode"), by.y = c("Aliquot", "Barcode"), all.x = T)
cellnumber_percluster_df <- barcode_merged_df %>%
  filter(is.na(predicted_doublet) | (!is.na(predicted_doublet) & predicted_doublet == F)) %>%
  mutate(id_cluster_uniq = Cluster_Name) %>%
  select(id_cluster_uniq) %>%
  table() %>%
  as.data.frame() %>%
  rename(id_cluster_uniq = '.') %>%
  mutate(easy_id = str_split_fixed(string = id_cluster_uniq, pattern = "_C", n = 2)[,1])
nrow(cellnumber_percluster_df[cellnumber_percluster_df$Freq >= 50,])
nrow(cellnumber_percluster_df[cellnumber_percluster_df$Freq > 50,])

clusters_selected_df <- cellnumber_percluster_df %>%
  filter(easy_id != "C3L-00359-T1") %>%
  filter(Freq > 50)
nrow(clusters_selected_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellNumberPerTumorManualCluster.", run_id, ".tsv")
write.table(x = cellnumber_percluster_df, file = file2write, quote = F, sep = "\t", row.names = F)
