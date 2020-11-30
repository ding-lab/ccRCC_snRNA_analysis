# Yige Wu @WashU Nov 2020

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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_each_atac_tumor_vs_pt_on_katmai/20201130.v2/findallmarkers_wilcox_each_snatac_tumor_vs_pt.20201130.v2.tsv")
## input samples to intersect
easyids2interesect <- c("C3L-00088-T1", "C3L-00088-T2", "C3L-00917-T1")

# summarize genes by occurance ---------------------------------------------
deg_filtered_df <- deg_df %>%
  filter(easyid_tumor %in% easyids2interesect) %>%
  filter(p_val_adj < 0.05)
deg_wide_df <- dcast(data = deg_filtered_df, formula = genesymbol_deg~easyid_tumor, value.var = "avg_logFC")
isallup_vec <- (deg_wide_df$`C3L-00088-T1`>0 & deg_wide_df$`C3L-00088-T2`>0 & deg_wide_df$`C3L-00917-T1`>0)
isalldown_vec <- (deg_wide_df$`C3L-00088-T1`<0 & deg_wide_df$`C3L-00088-T2`<0 & deg_wide_df$`C3L-00917-T1`<0)
deg_wide_df$direction_shared <- ifelse(isallup_vec, "up",
                                       ifelse(isalldown_vec, "down", NA))
deg_wide_df$mean_avg_logFC <- rowMeans(x = deg_wide_df[,easyids2interesect], na.rm = T)
deg_wide_df <- deg_wide_df %>%
  arrange(desc(direction_shared), desc(mean_avg_logFC))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs_nonmutant_snatac_tumors_shared", ".tsv")
write.table(x = deg_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)

