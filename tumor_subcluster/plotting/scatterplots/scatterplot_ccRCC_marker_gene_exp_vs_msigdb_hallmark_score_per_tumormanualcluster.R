# Yige Wu @WashU Apr 2021

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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input correaltion
spearman_result_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/calculate_correlation_for_ccRCC_markers_sct_data_exp_vs_pathwayscore_across_tumormanualclusters/20210712.v2/ccRCC_marker_expression_vs_pathway_score_bytumorcluster.correlation.tsv")
## input cell type fraction
pathscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_top_geneset_scores/20210419.v1/MSigDB.Hallmark.tsv")
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/avgeexp_tumorcells_sct_data_by_manualcluster_rm_doublets_on_katmai/20210413.v1/AverageExpression_ByManualTumorSubcluster.20210413.v1.tsv", data.table = F)
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv", data.table = F)

# format expression data --------------------------------------------------
spearman_result_plot_df <- spearman_result_df %>%
  filter(fdr < 0.05)
genes_process <- unique(spearman_result_plot_df$gene)
exp_filtered_long_df <- avgexp_df %>%
  filter(V1 %in% genes_process) %>%
  melt() %>%
  mutate(id_bycluster_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  mutate(easyid_column = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(cluster_name = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,2])


# make data frame for plotting --------------------------------------------
i <- 1
## make colors
colors_sample <- Polychrome::palette36.colors(n = length(unique(exp_filtered_long_df$easyid_column)))
names(colors_sample) <- unique(exp_filtered_long_df$easyid_column)
for (i in 1:nrow(spearman_result_plot_df)) {
  gene_tmp <- spearman_result_plot_df$gene[i]
  scorename_tmp <- spearman_result_plot_df$scorename[i]
  plotdata_df <- merge(x = pathscores_df[, c(scorename_tmp, "cluster_name")],
                       y = exp_filtered_long_df %>%
                         filter(V1 == gene_tmp),
                       by.x = c("cluster_name"), by.y = c("id_bycluster_byaliquot"))
  plotdata_df[, "value.y"] <- plotdata_df[, scorename_tmp]
  plotdata_df[, "value.x"] <- plotdata_df[, "value"]
  
  p <- ggplot()
  p <- ggscatter(data = plotdata_df, x = "value.x", y = "value.y", color = "easyid_column",
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE # Add confidence interval
  )
  # Add correlation coefficient
  p <- p + stat_cor(method = "spearman", 
                    label.x = quantile(x = plotdata_df$value.x, probs = 0.95, na.rm = T), 
                    label.y = quantile(x = plotdata_df$value.y, probs = 0.95, na.rm = T))
  p <- p + xlab(paste0(gene_tmp, " expression by tumor cluster"))
  p <- p + ylab(paste0(scorename_tmp, " by tumor cluster"))
  p <- p + scale_color_manual(values = colors_sample)
  p <- p + theme(legend.position = "bottom")
  p <- p + guides(color = guide_legend(title = NULL, ncol = 4, label.theme = element_text(size = 8)))
  file2write <- paste0(dir_out, gene_tmp, "_", scorename_tmp,".png")
  png(file2write, width = 800, height = 1000, res = 150)
  print(p)
  dev.off()
}

