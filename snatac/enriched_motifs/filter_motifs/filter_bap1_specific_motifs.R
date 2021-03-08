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
## input samples to intersect
easyids2interesect <- c("C3L-00416-T2", "C3L-01313-T1", "C3N-01200-T1")
easyids2interesect <- c("C3L-00416-T2", "C3L-01313-T1", "C3N-01200-T1", "C3L-01287-T1")

# summarize genes by occurance ---------------------------------------------
dam_filtered_df <- dam_df %>%
  mutate(easyid_tumor = gsub(x = cell_t2, pattern = "Tumor_", replacement = "")) %>%
  filter(easyid_tumor %in% easyids2interesect) %>%
  filter(FDR < 0.05) %>%
  filter(mean_score2 > 0 & diff > 0 | mean_score1 > 0 & diff < 0)
dam_wide_df <- dcast(data = dam_filtered_df, formula = TF_Name~easyid_tumor, value.var = "diff")
isallup_vec <- (dam_wide_df$`C3L-00416-T2`>0 & dam_wide_df$`C3L-01313-T1`>0 & dam_wide_df$`C3N-01200-T1`>0)
isalldown_vec <- (dam_wide_df$`C3L-00416-T2`<0 & dam_wide_df$`C3L-01313-T1`<0 & dam_wide_df$`C3N-01200-T1`<0)
dam_wide_df$direction_shared <- ifelse(isallup_vec, "up",
                                       ifelse(isalldown_vec, "down", NA))
dam_wide_df$mean_diff <- rowMeans(x = dam_wide_df[,easyids2interesect], na.rm = T)
dam_wide_df <- dam_wide_df %>%
  arrange(desc(direction_shared), desc(mean_diff))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DAMs_bap1mutant_snatac_tumors", ".tsv")
write.table(x = dam_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)

