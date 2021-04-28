# Yige Wu @WashU Aug 2020

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
## input barcode-to-cell-cycle-phase
barcode2phase_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/plotting/plot_cellcycleassigned_in_selected_tumors/20210414.v1/TumorCellPhase.DoubletRemoved.tsv")
## input barcode-to-manual-cluster
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20201130.v1/Barcode2TumorSubclusterId.20201130.v1.tsv", data.table = F)

# merge -------------------------------------------------------------------
barcode_merged_df <- merge(x = barcode2phase_df, y = barcode2tumorsubcluster_df, by.x = c("orig.ident", "barcode"), by.y = c("orig.ident", "barcode"), all.x = T)
## count
count_byphase_bycluster_df <- barcode_merged_df %>%
  select(Cluster_Name, Phase) %>%
  table() %>%
  as.data.frame()
count_bycluster_df <- barcode_merged_df %>%
  select(Cluster_Name) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::rename(Cluster_Name = '.')
count_merged_df <- merge(x = count_byphase_bycluster_df, y = count_bycluster_df, by = c("Cluster_Name"), all.x = T, suffixes = c(".phase.cluster", ".cluster"))
count_merged_df <- count_merged_df %>%
  mutate(Frac.phase.cluster = (Freq.phase.cluster/Freq.cluster))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, 'Fraction_cells.ByPhase.ByManualCluster.tsv')
write.table(x = count_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
