# Yige Wu @WashU Aug 2021

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
deg_male_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_wilcox_male_ccRCC_vs_pt_on_katmai/20210809.v1/Male.ccRCC.wilcox.logfc.threshold0.25.min.pct0.1.min.diff.pct0.1.AssayRNA.tsv")
deg_female_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_wilcox_female_ccRCC_vs_pt_on_katmai/20210809.v1/Female.ccRCC.wilcox.logfc.threshold0.25.min.pct0.1.min.diff.pct0.1.AssayRNA.tsv")

# combine -----------------------------------------------------------------
deg_long_df <- rbind(deg_male_df %>%
                       mutate(Sex = "M"),
                     deg_female_df %>%
                       mutate(Sex = "F"))
deg_logn_sig_df <- deg_long_df %>%
  filter(p_val_adj < 0.05)
deg_summary_bysex_df <- dcast(data = deg_logn_sig_df, formula = genesymbol_deg ~ Sex, value.var = "avg_log2FC")
colnames(deg_summary_bysex_df) <- c("genesymbol.human", "avg_log2FC.F", "avg_log2FC.M")
deg_summary_bysex_df <- deg_summary_bysex_df %>%
  mutate(deg_cat_bysex = ifelse(is.na(avg_log2FC.F), ifelse(avg_log2FC.M > 0, "Male_Up", "Male_Down"), 
                                ifelse(is.na(avg_log2FC.M), ifelse(avg_log2FC.F > 0, "Female_Up", "Female_Down"), 
                                       ifelse(avg_log2FC.F > 0, ifelse(avg_log2FC.M > 0, "Both_Up", "Female_Up|Male_Down"),
                                ifelse(avg_log2FC.F < 0, ifelse(avg_log2FC.M < 0, "Both_Down", "Male_Up|Female_Down"), "Other")))))
table(deg_summary_bysex_df$deg_cat_bysex)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_vs_PT.snRNA_DEGs.BySex.", run_id, ".tsv")
write.table(x = deg_long_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ccRCC_vs_PT.snRNA_DEGs.BySex.Summary.", run_id, ".tsv")
write.table(x = deg_summary_bysex_df, file = file2write, quote = F, sep = "\t", row.names = F)
