# Yige Wu @WashU Jun 2021

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
deg_ind_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/summarize_tumor_vs_pt_DEGs/20210429.v1/Tumor_DEGs.20210429.v1.tsv")
deg_cnvcor_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/format_alltumor_vs_pt_CNV_corrected_degs/20210608.v1/LR.logfc.threshold0.1.min.pct0.1.min.diff.pct0.1.AssayRNA.tsv")

# overlap -----------------------------------------------------------------
## preprocess CNV corrected degs
deg_merged_df <- merge(x = deg_ind_df, 
                       y = deg_cnvcor_df, by = c("genesymbol_deg"), all = T)

deg_filtered1_df <- deg_merged_df %>%
  filter((Num_sig_up >= 15 & Num_down == 0) | (Num_sig_down >= 15 & Num_up == 0))
deg_filtered2_df <- deg_filtered1_df %>%
  filter(FDR.CNVcorrected < 0.05) %>%
  filter((Num_sig_up >= 15 & avg_log2FC.allTumorcellsvsPT > 0) | (Num_sig_down >= 15 & avg_log2FC.allTumorcellsvsPT < 0))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Tumor_vs_PT_DEGs.with.CNVcorrection.", run_id, ".tsv")
write.table(x = deg_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Consistent.Tumor_vs_PT_DEGs.CNVcorrected.", run_id, ".tsv")
write.table(x = deg_filtered2_df, file = file2write, quote = F, sep = "\t", row.names = F)
