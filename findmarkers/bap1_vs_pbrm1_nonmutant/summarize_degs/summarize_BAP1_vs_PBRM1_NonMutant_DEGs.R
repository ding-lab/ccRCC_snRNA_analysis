# Yige Wu @WashU Mar 2021

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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/unite_BAP1_DEGs/20210407.v1/BAP1_DEGs.20210407.v1.tsv")

# summarize genes by occurance ---------------------------------------------
deg_sig_long_df <- deg_df %>%
  filter(group1_mut_category == "BAP1 mutated") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(deg_category = paste0("BAP1", "_", ifelse(avg_logFC > 0, "Up", "Down")))
deg_wide_df <- dcast(data = deg_sig_long_df, formula = genesymbol_deg~deg_category)
## including all fold changes
deg_long_df <- deg_df %>%
  filter(group1_mut_category == "BAP1 mutated")
deg_wide_df2 <- dcast(data = deg_long_df, formula = genesymbol_deg~easyid_tumor, value.var = 'avg_logFC', na.rm = T)
easyids_tumor <- unique(deg_long_df$easyid_tumor)
deg_wide_df2$Num_up <- rowSums(deg_wide_df2[, easyids_tumor] > 0, na.rm = T)
deg_wide_df2$Num_down <- rowSums(deg_wide_df2[, easyids_tumor] < 0, na.rm = T)
## including only significant fold changes
deg_wide_df3 <- dcast(data = deg_sig_long_df, formula = genesymbol_deg~easyid_tumor, value.var = 'avg_logFC', na.rm = T)
deg_wide_df3$mean_avg_logFC <- rowMeans(deg_wide_df3[, easyids_tumor], na.rm = T)

## combine
deg_wide_df <- merge(x = deg_wide_df, y = deg_wide_df2, by = c("genesymbol_deg"), all.x = T)
deg_wide_df$mean_avg_logFC <- mapvalues(x = deg_wide_df$genesymbol_deg, from = deg_wide_df3$genesymbol_deg, to = as.vector(deg_wide_df3$mean_avg_logFC))
deg_wide_df$mean_avg_logFC <- as.numeric(deg_wide_df$mean_avg_logFC)
# filter BAP1-specific DEGs ----------------------------------------------
# cutoff_bap1_vs_others <- 0.5*length(unique(deg_df$easyid_tumor[deg_df$comparison == "BAP1_vs_PBRM1_Mutants_Tumorcells"]))
cutoff_bap1_vs_others <- 0.4*length(easyids_tumor)
deg_wide_df <- deg_wide_df %>%
  mutate(BAP1_vs_OtherTumor_snRNA = ifelse(BAP1_Up >= cutoff_bap1_vs_others, "Up",
                                           ifelse(BAP1_Down >= cutoff_bap1_vs_others, "Down", "Inconsistent")))
table(deg_wide_df$BAP1_vs_OtherTumor_snRNA)
deg_bap1_df <- deg_wide_df %>%
  filter(BAP1_vs_OtherTumor_snRNA != "Inconsistent") %>%
  arrange(desc(BAP1_Up), desc(BAP1_Down))
nrow(deg_bap1_df)
deg_bap1_df <- deg_bap1_df %>%
  filter(Num_up == 0 | Num_down == 0) %>%
  arrange(desc(BAP1_Up), desc(BAP1_Down))
nrow(deg_bap1_df)
table(deg_bap1_df$BAP1_vs_OtherTumor_snRNA)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_DEGs.Sig.", run_id, ".tsv")
write.table(x = deg_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_DEGs.Consistent", run_id, ".tsv")
write.table(x = deg_bap1_df, file = file2write, quote = F, sep = "\t", row.names = F)
