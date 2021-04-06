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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/summarize_degs/unite_PBRM1_BAP1_DEGs/20210326.v1/PBRM1_BAP1_DEGs.20210326.v1.tsv")

# summarize genes by occurance ---------------------------------------------
deg_filtered_df <- deg_df %>%
  filter(p_val_adj < 0.05) %>%
  mutate(deg_category = paste0(comparison, "_", ifelse(avg_logFC > 0, "Up", "Down")))
deg_wide_df <- dcast(data = deg_filtered_df, formula = genesymbol_deg~deg_category)
deg_wide_df2 <- dcast(data = deg_filtered_df, formula = genesymbol_deg~deg_category, value.var = 'avg_logFC', fun.aggregate = mean, na.rm = T)


# filter BAP1-specific DEGs ----------------------------------------------
# cutoff_bap1_vs_others <- 0.5*length(unique(deg_df$easyid_tumor[deg_df$comparison == "BAP1_vs_PBRM1_Mutants_Tumorcells"]))
cutoff_bap1_vs_others <- 0.4*length(unique(deg_df$easyid_tumor[deg_df$comparison == "BAP1_vs_PBRM1_Mutants_Tumorcells"]))
deg_wide_df <- deg_wide_df %>%
  mutate(BAP1_vs_PBRM1_Mutants_snRNA = ifelse(BAP1_vs_PBRM1_Mutants_Tumorcells_Up >= cutoff_bap1_vs_others, "Up",
                                              ifelse(BAP1_vs_PBRM1_Mutants_Tumorcells_Down >= cutoff_bap1_vs_others, "Down", "Inconsistent"))) %>%
  mutate(BAP1_vs_NonMutants_snRNA = ifelse(BAP1_vs_NonMutants_Tumorcells_Up >= cutoff_bap1_vs_others, "Up",
                                           ifelse(BAP1_vs_NonMutants_Tumorcells_Down >= cutoff_bap1_vs_others, "Down", "Inconsistent"))) %>%
  mutate(BAP1_Tumorcells_vs_PTcells_snRNA = ifelse(BAP1_Tumorcells_vs_PTcells_Up >= cutoff_bap1_vs_others, "Up",
                                           ifelse(BAP1_Tumorcells_vs_PTcells_Down >= cutoff_bap1_vs_others, "Down", "Inconsistent")))
table(deg_wide_df$BAP1_Tumorcells_vs_PTcells_snRNA)
table(deg_wide_df$BAP1_vs_NonMutants_snRNA)
table(deg_wide_df$BAP1_vs_PBRM1_Mutants_snRNA)

deg_bap1_df <- deg_wide_df %>%
  filter((BAP1_vs_PBRM1_Mutants_snRNA == "Up" & BAP1_vs_NonMutants_snRNA == "Up" & BAP1_Tumorcells_vs_PTcells_snRNA == "Up") | (BAP1_vs_PBRM1_Mutants_snRNA == "Down" & BAP1_vs_NonMutants_snRNA == "Down" & BAP1_Tumorcells_vs_PTcells_snRNA == "Down")) %>%
  arrange(desc(BAP1_vs_PBRM1_Mutants_Tumorcells_Up), desc(BAP1_vs_NonMutants_Tumorcells_Up), desc(BAP1_Tumorcells_vs_PTcells_Up))
nrow(deg_bap1_df)
table(deg_bap1_df$BAP1_Tumorcells_vs_PTcells_snRNA)

# filter PBRM1-specific DEGs ----------------------------------------------
# cutoff_PBRM1_vs_others <- 0.5*length(unique(deg_df$easyid_tumor[deg_df$comparison == "PBRM1_vs_BAP1_Mutants_Tumorcells"]))
cutoff_pbrm1_vs_others <- 0.4*length(unique(deg_df$easyid_tumor[deg_df$comparison == "PBRM1_vs_BAP1_Mutants_Tumorcells"]))

deg_wide_df <- deg_wide_df %>%
  mutate(PBRM1_vs_BAP1_Mutants_snRNA = ifelse(PBRM1_vs_BAP1_Mutants_Tumorcells_Up >= cutoff_pbrm1_vs_others, "Up",
                                              ifelse(PBRM1_vs_BAP1_Mutants_Tumorcells_Down >= cutoff_pbrm1_vs_others, "Down", "Inconsistent"))) %>%
  mutate(PBRM1_vs_NonMutants_snRNA = ifelse(PBRM1_vs_NonMutants_Tumorcells_Up >= cutoff_pbrm1_vs_others, "Up",
                                            ifelse(PBRM1_vs_NonMutants_Tumorcells_Down >= cutoff_pbrm1_vs_others, "Down", "Inconsistent"))) %>%
  mutate(PBRM1_Tumorcells_vs_PTcells_snRNA = ifelse(PBRM1_Tumorcells_vs_PTcells_Up >= cutoff_pbrm1_vs_others, "Up",
                                                    ifelse(PBRM1_Tumorcells_vs_PTcells_Down >= cutoff_pbrm1_vs_others, "Down", "Inconsistent")))
table(deg_wide_df$PBRM1_Tumorcells_vs_PTcells_snRNA)
table(deg_wide_df$PBRM1_vs_NonMutants_snRNA)
table(deg_wide_df$PBRM1_vs_BAP1_Mutants_snRNA)

deg_pbrm1_df <- deg_wide_df %>%
  filter((PBRM1_vs_BAP1_Mutants_snRNA == "Up" & PBRM1_vs_NonMutants_snRNA == "Up" & PBRM1_Tumorcells_vs_PTcells_snRNA == "Up") | (PBRM1_vs_BAP1_Mutants_snRNA == "Down" & PBRM1_vs_NonMutants_snRNA == "Down" & PBRM1_Tumorcells_vs_PTcells_snRNA == "Down")) %>%
  arrange(desc(PBRM1_vs_BAP1_Mutants_Tumorcells_Up), desc(PBRM1_vs_NonMutants_Tumorcells_Up), desc(PBRM1_Tumorcells_vs_PTcells_Up))
nrow(deg_pbrm1_df)
table(deg_pbrm1_df$PBRM1_Tumorcells_vs_PTcells_snRNA)


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_PBRM1_DEGs.Num_samples.", run_id, ".tsv")
write.table(x = deg_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_DEGs.", run_id, ".tsv")
write.table(x = deg_bap1_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "PBRM1_DEGs.", run_id, ".tsv")
write.table(x = deg_pbrm1_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_PBRM1_DEGs.Mean_avg_logFC.", run_id, ".tsv")
write.table(x = deg_wide_df2, file = file2write, quote = F, sep = "\t", row.names = F)
