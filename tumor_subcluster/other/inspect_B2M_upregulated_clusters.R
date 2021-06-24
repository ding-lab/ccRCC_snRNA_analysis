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
## input DEGs
deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/unite_degs/unite_degs_for_tumor_manualcluster/20210413.v1/TumorManualCluster.DEGs.Wilcox.Minpct0.1.Logfc0.min.diff.pct0.1.tsv")
## input cell-cell interactions
interactions_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/filter_interactions/filter_cellphonedb_out/20201012.v1/cell.phone.res.total.run20200818.filtered.txt")
interactions_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/annotate_interactions/annotate_cellphone_out/20201012.v1/cell.phone.res.total.run20200818.filtered.formatted.txt")

# filter ------------------------------------------------------------------
deg_filtered_df <- deg_all_df %>%
  filter(gene == "B2M") %>%
  filter(p_val_adj < 0.05) %>%
  filter(avg_logFC > 0)

interactions_b2m_df <- interactions_df %>%
  filter(gene.source == "B2M")

