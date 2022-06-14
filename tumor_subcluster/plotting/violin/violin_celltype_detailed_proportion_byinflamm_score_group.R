# Yige Wu @WashU June 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "data.table",
  "stringr",
  "ggpubr",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input data --------------------------------------------------------------
## input cell type fraction
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")
## input gene sets to test
genesets_test_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/pathway/count_ora_sig_genesets_in_up_degs_across_samples/20220606.v1/Count_gene_set_in_up_tumorcluster_degs.20220606.v1.tsv")
## input cell type proportion
celltype_frac_long_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/count_fraction/count_celltypedetailed_fraction_per_sample/20220607.v1/CellGroupBarcodes_Number_and_Fraction_per_Sample20220607.v1.tsv")

# pre-process -------------------------------------------------------------
genesets_test <- genesets_test_df$Description
scoregroups_test <- paste0(gsub(x = genesets_test_df$Description, pattern = "HALLMARK_", replacement = ""), "_Score")
cellgroups_test <- unique(celltype_frac_long_df$Cell_group); cellgroups_test <- cellgroups_test[!(cellgroups_test %in% c("Unknown", "Mixed myeloid/lymphoid", 
                                                                                                                         "CD4/CD8 proliferating"))]; cellgroups_test
celltype_frac_wide_df <- dcast(data = celltype_frac_long_df, formula = Aliquot_WU~Cell_group, value.var = "Frac_CellGroupBarcodes_ByAliquot")
celltype_frac_wide_df[is.na(celltype_frac_wide_df)] <- 0

## make plot data
scores_bycluster_df <- scores_df[, c("cluster_name", "INFLAMMATORY_RESPONSE_Score")]
colnames(scores_bycluster_df) <- c("cluster_name", "score.bycluster")
scores_bycluster_df <- scores_bycluster_df %>%
  mutate(sampleid = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
  mutate(sampleid = gsub(x = sampleid, pattern = "\\.", replacement = "-"))
sampleids_top_inflamm <- scores_bycluster_df$sampleid[scores_bycluster_df$score.bycluster >= quantile(x = scores_bycluster_df$score.bycluster, probs = 0.9)]
sampleids_bottom_inflamm <- scores_bycluster_df$sampleid[scores_bycluster_df$score.bycluster <= quantile(x = scores_bycluster_df$score.bycluster, probs = 0.1)]
plot_data_comp_df <- celltype_frac_wide_df %>%
  mutate(sample_group = ifelse(Aliquot_WU %in% sampleids_top_inflamm, "top",
                               ifelse(Aliquot_WU %in% sampleids_bottom_inflamm, "bottom", "middle")))

# test --------------------------------------------------------------------
my_comparisons <- list(c("top", "bottom"), c("top", "middle"))
for (cellgroup_tmp in cellgroups_test) {
  plotdata_df <- plot_data_comp_df[, c("Aliquot_WU", "sample_group", cellgroup_tmp)]
  colnames(plotdata_df) <- c("Aliquot_WU",  "x_plot", "y_plot")
  ## plot
  p <- ggviolin(data = plotdata_df, x = "x_plot", y = "y_plot",
                 add = "dotplot"
  )
  p <- p + stat_compare_means(method = "t.test", comparisons = my_comparisons)
  # p <- p + ggtitle(label = paste0(protein_tmp, "~", treatment_tmp, " 1-month response"))
  p <- p + theme_classic(base_size = 17)
  p <- p + xlab(paste0("inflammation score"))
  p <- p + ylab(paste0("% ", cellgroup_tmp))
  # p <- p + xlim(c(min(plotdata_tmp_df$x_plot)-0.11, max(plotdata_tmp_df$x_plot)+0.11))
  # p <- p + ylim(c(min(plotdata_tmp_df$y_plot)-0.05, max(plotdata_tmp_df$y_plot)+0.1))
  p <- p + theme(axis.text = element_text(color = "black"))
  # p
  ## write output
  file2write <- paste0(dir_out, gsub(x = cellgroup_tmp, pattern = " |\\/", replacement = "_"), ".png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()
}


