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
deg_filtered_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/filter_deg/filter_each_tumor_vs_pt_degs_bysnatactumorgroup_shared/20201202.v1/Top_DEGs_Tumorcells_vs_PT_ByTumorGroup.20201202.v1.tsv")
deg_avg_loFC_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_deg_by_snatactumor_avg_logFC/20201204.v1/DEGs_all_snatac_tumors_avg_logFC.20201204.v1.tsv")

# merge -------------------------------------------------------------------
deg_merged_df <- merge(x = deg_filtered_df, y = deg_avg_loFC_df, by = c("genesymbol_deg"), all.x = T)
deg_merged_df <- deg_merged_df %>%
  arrange(category_byshared_label, desc(mean_avg_logFC.bap1mutant), desc(mean_avg_logFC.pbrm1mutant), desc(mean_avg_logFC.nonmutant))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Filtered_degs_snatac_tumors_avg_logFC.", run_id, ".tsv")
write.table(x = deg_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
