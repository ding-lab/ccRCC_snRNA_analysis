# Yige Wu @WashU Dec 2020

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

# input dependencies ------------------------------------------------------
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210305.v1/meta_data.20210305.v1.tsv")
## input mutation info
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/merge_bulk_sn_profiles/20210504.v1/bulk_sn_omics_profile.20210504.v1.tsv")
## input cell type fraction
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")
## input gene sets to test
genesets_test_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/pathway/count_ora_sig_genesets_in_up_degs_across_samples/20220606.v1/Count_gene_set_in_up_tumorcluster_degs.20220606.v1.tsv")

# set plotting parameters -------------------------------------------------
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.25, seed = 2)
genes_test <- c("BAP1", "KDM5C", "PBRM1", "SETD2")

# preprocess --------------------------------------------------------------
scores_df <- scores_df %>%
  mutate(sample_id = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
  mutate(sample_id = gsub(x = sample_id, pattern = "\\.", replacement = "-")) %>%
  filter(sample_id != "C3L-00359-T1")
scoregroups_test <- paste0(gsub(x = genesets_test_df$Description, pattern = "HALLMARK_", replacement = ""), "_Score")

# process -----------------------------------------------------------------
result_df <- NULL
for (gene_mut in genes_test) {
  mut_tmp_df <- mut_df[, c("Aliquot_snRNA_WU", paste0("Mut.", gene_mut))]
  colnames(mut_tmp_df) <- c("Aliquot_snRNA_WU", "mutation_class")
  
  wilcox_pvalues_vec <- NULL
  ttest_pvalue_vec <- NULL
  for (scoregroup_tmp in scoregroups_test) {
    ## make plot data
    scores_tmp_df <- scores_df[, c("cluster_name", "sample_id", scoregroup_tmp)]
    colnames(scores_tmp_df) <- c("cluster_name", "sample_id", "score.bycluster")
    scores_bysample_df <- scores_tmp_df %>%
      group_by(sample_id) %>%
      summarise(score.bysample = max(score.bycluster))
    
    plot_data_df <- merge(x = scores_bysample_df,
                          y = mut_tmp_df,
                          by.x = c("sample_id"), by.y = c("Aliquot_snRNA_WU"), all.x = T)
    plot_data_df <- plot_data_df %>%
      filter(!is.na(mutation_class) & !is.na(score.bysample)) %>%
      mutate(sample_group = ifelse(mutation_class == "None", "WT", "Mutated")) %>%
      mutate(y_plot = score.bysample)
    plot_data_df$x_plot <- factor(x = plot_data_df$sample_group, levels = c("WT", "Mutated"))
    # ## plot
    # p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
    # p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 0.6)
    # p = p + geom_boxplot(width=.1, outlier.shape = 23, outlier.fill = "red")
    # p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
    # p = p + stat_compare_means(data = plot_data_df, 
    #                            mapping = aes(x = x_plot, y = y_plot), ref.group = "WT",
    #                            symnum.args = symnum.args,
    #                            hide.ns = F, method = "wilcox.test")
    # p <- p + theme_classic()
    # p <- p + theme(legend.position = "none")
    # p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 9.5, colour = "black"))
    # p <- p + ylab(label = scoregroup_tmp)
    # file2write <- paste0(dir_out, gene_mut, "_", gsub(x = scoregroup_tmp, pattern = "\\/", replacement = "_"), ".png")
    # png(file2write, width = 600, height = 600, res = 150)
    # print(p)
    # dev.off()
    
    ## test
    ## wilcox test
    test_result <- wilcox.test(plot_data_df$y_plot[plot_data_df$x_plot == "Mutated"], plot_data_df$y_plot[plot_data_df$x_plot == "WT"])
    wilcox_pvalues_vec <- c(wilcox_pvalues_vec, test_result$p.value)
    ## t test
    test_result <- t.test(plot_data_df$y_plot[plot_data_df$x_plot == "Mutated"], plot_data_df$y_plot[plot_data_df$x_plot == "WT"])
    ttest_pvalue_vec <- c(ttest_pvalue_vec, test_result$p.value)
  }
  result_tmp_df <- data.frame(genesymbol_mut = gene_mut, name_score = scoregroups_test, 
                              pvalue.wilcox = wilcox_pvalues_vec, fdr.wilcox.bymutgene = p.adjust(p = wilcox_pvalues_vec, method = "fdr"),
                              pvalue.ttest = ttest_pvalue_vec, fdr.ttest.bymutgene = p.adjust(p = ttest_pvalue_vec, method = "fdr"))
  result_df <- rbind(result_df, result_tmp_df)
}
result_df$fdr.wilcox = p.adjust(p = result_df$pvalue.wilcox, method = "fdr")
result_df$fdr.ttest = p.adjust(p = result_df$pvalue.ttest, method = "fdr")

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "wilcox_ttest_results.tsv")
write.table(x = result_df, file = file2write, quote = F, sep = "\t", row.names = F)


# plot --------------------------------------------------------------------
# process -----------------------------------------------------------------
for (gene_mut in "BAP1") {
  mut_tmp_df <- mut_df[, c("Aliquot_snRNA_WU", paste0("Mut.", gene_mut))]
  colnames(mut_tmp_df) <- c("Aliquot_snRNA_WU", "mutation_class")
  
  wilcox_pvalues_vec <- NULL
  ttest_pvalue_vec <- NULL
  for (scoregroup_tmp in c("PROTEIN_SECRETION_Score", "MTORC1_SIGNALING_Score")) {
    ## make plot data
    scores_tmp_df <- scores_df[, c("cluster_name", "sample_id", scoregroup_tmp)]
    colnames(scores_tmp_df) <- c("cluster_name", "sample_id", "score.bycluster")
    scores_bysample_df <- scores_tmp_df %>%
      group_by(sample_id) %>%
      summarise(score.bysample = max(score.bycluster))
    
    plot_data_df <- merge(x = scores_bysample_df,
                          y = mut_tmp_df,
                          by.x = c("sample_id"), by.y = c("Aliquot_snRNA_WU"), all.x = T)
    plot_data_df <- plot_data_df %>%
      filter(!is.na(mutation_class) & !is.na(score.bysample)) %>%
      mutate(sample_group = ifelse(mutation_class == "None", "WT", "Mutated")) %>%
      mutate(y_plot = score.bysample)
    plot_data_df$x_plot <- factor(x = plot_data_df$sample_group, levels = c("WT", "Mutated"))
    ## plot
    p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
    p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 0.6)
    p = p + geom_boxplot(width=.1, outlier.shape = 23, outlier.fill = "red")
    p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
    p = p + stat_compare_means(data = plot_data_df,
                               mapping = aes(x = x_plot, y = y_plot), ref.group = "WT",
                               symnum.args = symnum.args,
                               hide.ns = F, method = "t.test")
    p <- p + theme_classic()
    p <- p + theme(legend.position = "none")
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 9.5, colour = "black"))
    p <- p + ylab(label = scoregroup_tmp)
    file2write <- paste0(dir_out, gene_mut, "_", gsub(x = scoregroup_tmp, pattern = "\\/", replacement = "_"), ".png")
    png(file2write, width = 600, height = 600, res = 150)
    print(p)
    dev.off()
  }
}


