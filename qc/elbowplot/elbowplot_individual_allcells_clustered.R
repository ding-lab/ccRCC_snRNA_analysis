# Yige Wu @WashU Dec 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
library(ggplot2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
pct_df <- fread(data.table = F, input = "./Resources/Analysis_Results/qc/elbowplot/elbowplot_individual_allcells_clustered_katmai/20220221.v1/pct_standard_variances_by_PC.20220221.v1.tsv")

# plot by sample ----------------------------------------------------------
easyid_tmp <- "C3L-01313-T1"
for (easyid_tmp in unique(pct_df$easy_id)) {
  plotdata_df <- pct_df %>%
    filter(easy_id == easyid_tmp)
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, mapping = aes(x = rank_pc, y = pct), size = 0.8)
  p <- p + theme_classic()
  p <- p + xlab("PC (ranked)") + ylab("Standard Deviation (%)")
  p <- p + ggtitle(label = paste0("Sample ", easyid_tmp))
  file2write <- paste0(dir_out, easyid_tmp, ".png")
  png(file2write, width = 500, height = 450, res = 150)
  print(p)
  dev.off()
}
