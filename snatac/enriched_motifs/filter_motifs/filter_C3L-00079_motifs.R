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
## input motif result
### mean score 1 is for PT, mean score 2 is for the tumor cells
dam_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Enriched_Motifs/Tumor_vs_NormalPT/Score_difference.Tumor_Normal_comparison.20201130.tsv")

# summarize genes by occurance ---------------------------------------------
dam_filtered_df <- dam_df %>%
  mutate(easyid_tumor = gsub(x = cell_t2, pattern = "Tumor_", replacement = "")) %>%
  filter(easyid_tumor %in% c("C3L-00079-T1")) %>%
  # filter(mean_score2 > 0 & diff > 0)  %>%
  filter(FDR < 0.05) 
dam_wide_df <- dcast(data = dam_filtered_df, formula = TF_Name~easyid_tumor, value.var = "diff")
countup_vec <- rowSums(x = (dam_wide_df[,-1] > 0), na.rm = T)
table(countup_vec)
dam_wide_filtered_df <- dam_wide_df[countup_vec >= 8,]
dam_wide_filtered_df <- dam_wide_filtered_df %>%
  dplyr::filter(!(TF_Name %in% c("NFIA", "NFIC::TLX1", "NFIC(var.2)", "NFIX(var.2)"))) %>% ## vecause the NFIC was already chosen and NFIA looks absent from non-mutant tumors
  dplyr::filter(!(TF_Name %in% c("HOXA2", "HOXD3", "GSX2", "EMX1"))) ## these are removed because of they are all similar but HOXB2 was reported to be up-regulated in cancer
  
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DAMs_all_snatac_tumors_shared", ".tsv")
write.table(x = dam_wide_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "DAMs_diff_acroos_all_snatac_tumors", ".tsv")
write.table(x = dam_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)

