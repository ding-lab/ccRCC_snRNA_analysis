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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/summarize_degs/summarize_PBRM1_BAP1_DEGs/20210405.v1/BAP1_PBRM1_DEGs.Num_samples.20210405.v1.tsv")

for (direction_deg in c("Up", "Down")) {
  # make plot data ----------------------------------------------------------
  ## filter
  plot_data_df <- deg_df %>%
    select(BAP1_Tumorcells_vs_PTcells_snRNA, BAP1_vs_PBRM1_Mutants_snRNA, BAP1_vs_NonMutants_snRNA)
  plot_data_df <- as.data.frame(plot_data_df == direction_deg)
  plot_data_df <- plot_data_df %>%
    filter(BAP1_Tumorcells_vs_PTcells_snRNA | BAP1_vs_PBRM1_Mutants_snRNA | BAP1_vs_NonMutants_snRNA)
  # venn --------------------------------------------------------------------
  plot_data_euler <- euler(plot_data_df)
  file2write <- paste0(dir_out, "BAP1", ".", direction_deg, ".png")
  png(file2write, height = 1000, width = 800, res = 150)
  p <- plot(plot_data_euler, labels = F, quantities = list(fontsize = 10), legend = T)
  grid.draw(p)
  dev.off()
}

