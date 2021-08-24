# Yige Wu @WashU Apr 2021

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
# deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_pt_DEGs/20210429.v1/Tumor_vs_PT_DEGs.20210429.v1.tsv")
deg_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/Tumorcells_vs_PTcells/20210823.v1/Each31RCC_Tumorcells_vs_CombinedPTcells.20210823.v1.tsv")

# summarize genes by occurance ---------------------------------------------
deg_sig_long_df <- deg_df %>%
  filter(easyid_tumor != "C3L-00359-T1") %>%
  filter(p_val_adj < 0.05) %>%
  mutate(deg_category = paste0("Num_sig_", ifelse(avg_log2FC > 0, "up", "down")))
deg_wide_df <- reshape2::dcast(data = deg_sig_long_df, formula = genesymbol_deg~deg_category)
## including all fold changes
deg_long_df <- deg_df %>%
  filter(easyid_tumor != "C3L-00359-T1") 
deg_wide_df2 <- reshape2::dcast(data = deg_long_df, formula = genesymbol_deg~easyid_tumor, value.var = 'avg_log2FC', na.rm = T)
easyids_tumor <- unique(deg_long_df$easyid_tumor)
deg_wide_df2$Num_up <- rowSums(deg_wide_df2[, easyids_tumor] > 0, na.rm = T)
deg_wide_df2$Num_down <- rowSums(deg_wide_df2[, easyids_tumor] < 0, na.rm = T)
deg_wide_df2$Num_datapoints <- rowSums(!is.na(deg_wide_df2[, easyids_tumor]), na.rm = T)
## including only significant fold changes
deg_wide_df3 <- reshape2::dcast(data = deg_sig_long_df, formula = genesymbol_deg~easyid_tumor, value.var = 'avg_log2FC', na.rm = T)
deg_wide_df3$mean_avg_log2FC <- rowMeans(deg_wide_df3[, easyids_tumor], na.rm = T)

## combine
deg_wide_df <- merge(x = deg_wide_df, y = deg_wide_df2, by = c("genesymbol_deg"), all.x = T)
deg_wide_df$mean_avg_log2FC <- mapvalues(x = deg_wide_df$genesymbol_deg, from = deg_wide_df3$genesymbol_deg, to = as.vector(deg_wide_df3$mean_avg_log2FC))
deg_wide_df$mean_avg_log2FC <- as.numeric(deg_wide_df$mean_avg_log2FC)
deg_wide_df <- deg_wide_df[, c("genesymbol_deg", "Num_sig_up", "Num_sig_down", "mean_avg_log2FC", "Num_datapoints", "Num_up", "Num_down", easyids_tumor)]

# filter DEGs ----------------------------------------------
cutoff_datapoints <- length(easyids_tumor)*0.5
deg_enough_datapoint_df <- deg_wide_df %>%
  filter(Num_datapoints > cutoff_datapoints) %>%
  mutate(Tumor_vs_PT = ifelse(Num_sig_up > cutoff_datapoints & Num_sig_down == 0, "Up",
                              ifelse(Num_sig_down > cutoff_datapoints & Num_sig_up == 0, "Down", "Inconsistent")))
table(deg_enough_datapoint_df$Tumor_vs_PT)
deg_consistent_df <-  deg_enough_datapoint_df %>%
  filter(Tumor_vs_PT != "Inconsistent")
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Tumor_DEGs.", run_id, ".tsv")
write.table(x = deg_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Tumor_DEGs.EnoughDataPoints.", run_id, ".tsv")
write.table(x = deg_enough_datapoint_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Tumor_DEGs.EnoughDataPoints.Consistent.", run_id, ".tsv")
write.table(x = deg_consistent_df, file = file2write, quote = F, sep = "\t", row.names = F)
