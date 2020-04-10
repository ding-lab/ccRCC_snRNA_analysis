# Yige Wu @WashU Apr 2020
## running on local
## for scaling average expression across cells

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression calculated (RNA)
avgexp_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/averageexpression/averageexpression_tumor_cells_by_manual_subcluster/20200325.v1/averageexpression_tumor_cells_by_manual_subcluster.20200325.v1.tsv", data.table = F)

# scale -------------------------------------------------------------------
avgexp_mat <- avgexp_df %>%
  select(-V1)
rownames(avgexp_mat) <- avgexp_df$V1
avgexp_mat <- as.matrix(avgexp_mat)
avgexp_mat_t <- t(avgexp_mat)
scaled_avgexp_mat_t <- scale(x = avgexp_mat_t, scale = T, center = T)
scaled_avgexp_mat <- t(scaled_avgexp_mat_t)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "averageexpression_tumor_cells_by_manual_subcluster.", "scaled.", run_id, ".tsv")
write.table(x = scaled_avgexp_mat, file = file2write, quote = F, sep = "\t", row.names = T)

