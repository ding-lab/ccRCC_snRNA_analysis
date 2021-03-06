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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_macrophage_deg_by_snatactumor_avg_logFC/20201215.v1/DEGs,macrophages.tumor_vs_nat.avg_logFC.20201215.v1.tsv")

# summarize genes by occurance ---------------------------------------------
deg_filtered_df <- deg_df %>%
  filter(count_up_tumors >=6) %>%
  arrange(desc(count_up_tumors), desc(count_down_tumors))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Up.DEGs.macrophages.tumor_vs_nat.avg_logFC.shared.", run_id, ".tsv")
write.table(x = deg_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

