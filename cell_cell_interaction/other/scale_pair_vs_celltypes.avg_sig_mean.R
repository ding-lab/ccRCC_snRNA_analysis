# Yige Wu @WashU Sep 2020
## reference of the cellphonedb output: https://www.cellphonedb.org/documentation

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
## input cellphonedb output
cellphone_long_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/annotate_interactions/annotate_cellphone_out/20200923.v1/cell.phone.res.total.run20200818.filtered.formatted.txt")

# count significant cases, average mean values across samples -----------------------------
table(cellphone_long_df$variable) %>% length()
cellphone_long_sum_df <- cellphone_long_df %>%
  group_by(interacting_pair, variable) %>%
  summarize(avg_sig_mean = mean(x = sig_mean, na.rm = T))
cellphone_wide_df <- dcast(data = cellphone_long_df, formula = pair_cell.types ~ Easy_id, value.var = "value")
cellphone_wide_mat <- as.matrix(cellphone_wide_df[,-1])
rownames(cellphone_wide_mat) <- cellphone_wide_df$pair_cell.types
cellphone_wide_mat[1:5, 1:10]
cellphone_scaled_mat <- scale(x = cellphone_wide_mat, center = T, scale = T)
cellphone_scaled_mat[1:5, 1:10]

summary(cellphone_scaled_mat[,1])

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cellphonedb.pair_cell.types_vs_sample.sig_mean", ".tsv")
write.table(x = cellphone_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "cellphonedb.pair_cell.types_vs_sample.sig_mean.scaled_by_sample.", ".tsv")
write.table(x = cellphone_scaled_mat, file = file2write, quote = F, sep = "\t", row.names = T)
