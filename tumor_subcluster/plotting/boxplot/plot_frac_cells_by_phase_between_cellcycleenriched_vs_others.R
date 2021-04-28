# Yige Wu @WashU Aug 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggpubr)

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input fraction of cells
frac_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/count/count_cellnumber_by_cellcyclephase_for_tumor_manualcluster/20210427.v1/Fraction_cells.ByPhase.ByManualCluster.tsv")
## input cluster-to-module-assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210421.v3/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

# set parameters ----------------------------------------------------------
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.25, seed = 2)

# make plot data ----------------------------------------------------------
plotdata_df <- frac_df %>%
  filter(Freq.cluster >= 50)
enrich_df <- enrich_df %>%
  mutate(cluster_name2 = gsub(x = cluster_name, pattern = "\\.", replacement = "-"))
plotdata_df <- plotdata_df %>%
  mutate(x_plot = ifelse(Cluster_Name %in% enrich_df$cluster_name2[enrich_df$Cell_cycle], "CellCycle\nEnriched", "Others")) %>%
  mutate(y_plot = Frac.phase.cluster)

# plot --------------------------------------------------------------------
for (phase_tmp in unique(plotdata_df$Phase)) {
  plotdata_tmp_df <- subset(x = plotdata_df, Phase == phase_tmp)
  ## plot
  p = ggplot(data = plotdata_tmp_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 0.6)
  p = p + geom_boxplot(width=.1, outlier.shape = 23, outlier.fill = "red")
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
  p = p + stat_compare_means(data = plotdata_tmp_df, 
                             mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..),  ref.group = "Others",
                             symnum.args = symnum.args,
                             hide.ns = F, method = "wilcox.test")
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 9.5, colour = "black"))
  p <- p + ylab(label = paste0("% of tumor cells in ", phase_tmp, " phase"))
  file2write <- paste0(dir_out, phase_tmp, ".png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()
}

