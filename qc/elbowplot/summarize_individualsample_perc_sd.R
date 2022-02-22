# Yige Wu @WashU Feb 2022
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
# ## set run id
# version_tmp <- 1
# run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
# ## set output directory
# dir_out <- paste0(makeOutDir(), run_id, "/")
# dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the % of standard deviation contributed by each PC
pct_long_df <- fread(data.table = F, input = "./Resources/Analysis_Results/qc/elbowplot/elbowplot_individual_allcells_clustered_katmai/20220221.v1/pct_standard_variances_by_PC.20220221.v1.tsv")

# reshape -----------------------------------------------------------------
pct_pc30_df <- pct_long_df %>%
  filter(rank_pc == 30)
summary(pct_pc30_df$cumu_pct)

pct_wide_df <- dcast(data = pct_long_df, )

