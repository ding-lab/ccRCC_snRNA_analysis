# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
# library(org.Hs.eg.db)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input degs
deg_79_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_C3L-00079_tumorlike_vs_tumorcells/20200909.v1/findmarkers_wilcox_tumorlikecells_vs_tumorcells.logfcthreshold0.693147180559945.minpct0.1.mindiffpct0.1.tsv")
deg_79_df$Case <- "C3L-00079"
deg_12_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_C3N-01200_vim_high_vs_low/20200910.v1/findmarkers_wilcox_vimhigh_vs_low.logfcthreshold0.693147180559945.minpct0.1.mindiffpct0.1.tsv")
deg_12_df <- deg_12_df %>%
  filter(clusterid_group1 == 6)
deg_12_df$Case <- "C3N-01200"

# merge -------------------------------------------------------------------
deg_merged_df <- merge(deg_79_df, deg_12_df, by = c("row_name"), all = T, suffixes = c(".79", ".12"))
deg_merged_df <- deg_merged_df %>%
  dplyr::rename(gene_symbol = row_name)
deg_shared_df <- deg_merged_df %>%
  dplyr::filter(!is.na(p_val_adj.79) & !is.na(p_val_adj.12)) %>%
  mutate(Direction = ifelse(p_val_adj.79 < 0.05 & p_val_adj.12 < 0.05, 
                            ifelse(avg_logFC.79 > 0 & avg_logFC.12 > 0,
                                   "Up-regulated",
                                   ifelse(avg_logFC.79 < 0 & avg_logFC.12 < 0,
                                          "Down-regulated",
                                          "inconsistent direction")),
                            "Not significant"))


# write output ------------------------------------------------------------


