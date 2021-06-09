# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
dam_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Enriched_Motifs/Score_difference.Tumor_Normal_comparison.20210509.tsv")

# make data for plotting --------------------------------------------------
dam_sig_df <- dam_df %>%
  mutate(log10_pvalue = -log10(pvalue)) %>%
  mutate(log10_pvalue_capped = ifelse(is.infinite(log10_pvalue), 150, log10_pvalue)) %>%
  filter(FDR < 0.05)

dam_sum_df <- dam_sig_df %>%
  group_by(TF_Name) %>%
  summarise(Num_sig_up = length(which(diff > 0)), Num_sig_down = length(which(diff < 0)),
            avg_log10_pvalue = mean(log10_pvalue_capped), avg_diff = mean(diff)) %>%
  mutate(foldchange_type = ifelse(Num_sig_down == 0, "consistently higher in ccRCC",
                                  ifelse(Num_sig_up == 0, "consistently lower in ccRCC", "mixed fold change directions")))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Summarized.ccRCC_vs_PT.Motifs.", run_id, ".tsv")
write.table(x = dam_sum_df, file = file2write, quote = F, sep = "\t", row.names = F)
