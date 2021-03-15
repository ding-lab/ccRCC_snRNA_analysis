# Yige Wu @WashU Mar 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
if (!requireNamespace("eulerr", quietly = TRUE))
  install.packages("eulerr")
library(eulerr)
library(grid)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input DEGs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/run_deg_analysis/run_PBRM1_BAP1_deg_on_cptac_ccRCC_discovery_cases/20210312.v2/PBRM1_BAP1_DEGs.glmQLFTest.outputtables.tsv")

# make plot data ----------------------------------------------------------
## filter
deg_filtered_df <- deg_df %>%
  filter(logFC > 0) %>%
  filter(FDR < 0.05) %>%
  filter(grepl(pattern = "vs Non-mutants", x = comparison)) %>%
  filter(!(gene_ensembl_id %in% c("__alignment_not_unique", "__ambiguous", "__no_feature")))
## make data frame
logFC_wide_df <- dcast(data = deg_filtered_df, formula = gene_ensembl_id ~ comparison, value.var = "logFC")
plot_data_df <- as.data.frame(!is.na(logFC_wide_df[, 2:3]))

# venn --------------------------------------------------------------------
plot_data_euler <- euler(plot_data_df)
file2write <- paste0(dir_out, "Up_DEGs.png")
png(file2write, height = 1000, width = 800, res = 150)
p <- plot(plot_data_euler, labels = F, quantities = list(fontsize = 10), legend = T)
grid.draw(p)
dev.off()
