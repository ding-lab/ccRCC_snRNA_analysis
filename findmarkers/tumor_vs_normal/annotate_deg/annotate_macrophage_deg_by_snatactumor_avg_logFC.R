# Yige Wu @WashU Dec 2020

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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_each_atac_tumor_vs_nat_macrophage_on_katmai/20201215.v1/findallmarkers_wilcox_each_snatac_tumor_vs_nat_macrophages.20201215.v1.tsv")

# summarize genes by occurance ---------------------------------------------
deg_filtered_df <- deg_df %>%
  filter(p_val_adj < 0.05)
deg_wide_df <- dcast(data = deg_filtered_df, formula = genesymbol_deg~easyid_tumor, value.var = "avg_logFC")
## summarize
easyids_tumor <- unique(deg_df$easyid_tumor)
deg_wide_df$count_up_tumors <- rowSums(x = (deg_wide_df[,easyids_tumor] > 0), na.rm = T)
deg_wide_df$count_down_tumors <- rowSums(x = (deg_wide_df[,easyids_tumor] < 0), na.rm = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs,macrophages.tumor_vs_nat.avg_logFC.", run_id, ".tsv")
write.table(x = deg_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)

