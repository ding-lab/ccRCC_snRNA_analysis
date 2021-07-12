# Yige Wu @WashU July 2021

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
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210702.v1/ccRCC_markers.Surface.20210702.v1.tsv")
## input druggable genes
druggenes_df <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/approved_target_ids_all.csv")

# filter ------------------------------------------------------------------
markers_merged_df <- merge(x = markers_df, y = druggenes_df, by.x = c("Gene"), by.y = c("Gene Name"), all.x = T)
markers_merged_df <- markers_merged_df %>%
  arrange(`Drug IDs`)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_markers.Surface.", run_id, ".tsv")
write.table(x = markers_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ccRCC_markers.Surface.", run_id, ".csv")
write.table(x = markers_merged_df, file = file2write, quote = F, sep = ",", row.names = F)
