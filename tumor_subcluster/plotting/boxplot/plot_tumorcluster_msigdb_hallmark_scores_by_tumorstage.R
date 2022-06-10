# Yige Wu @WashU Apr 2021

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

# input  ------------------------------------------------------
## input clinical data
clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_case_clinical_data/20201125.v1/snRNA_ccRCC_Clinicl_Table.20201125.v1.tsv")
## input cell type fraction
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")
## input gene sets to test
genesets_test_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/pathway/count_ora_sig_genesets_in_up_degs_across_samples/20220606.v1/Count_gene_set_in_up_tumorcluster_degs.20220606.v1.tsv")


# preprocess -------------------------------------------------
pos <- position_jitter(width = 0.2, seed = 1)
clinical_filtered_df <- clinical_df %>% 
  filter(Case != "C3L-00359") %>%
  select(Case, Tumor_Stage_Pathological) 
# scoregroups_test <- colnames(scores_df)
# scoregroups_test <- scoregroups_test[grepl(pattern = "Score", x = scoregroups_test)]
scoregroups_test <- c("UV_RESPONSE_DN_Score", "HYPOXIA_Score", "MITOTIC_SPINDLE_Score",
                      "INTERFERON_GAMMA_RESPONSE_Score", "COMPLEMENT_Score", "EPITHELIAL_MESENCHYMAL_TRANSITION_Score",
                      "E2F_TARGETS_Score", "KRAS_SIGNALING_UP_Score", "MTORC1_SIGNALING_Score",
                      "ALLOGRAFT_REJECTION_Score", "G2M_CHECKPOINT_Score", "INFLAMMATORY_RESPONSE_Score",
                      "MYC_TARGETS_V1_Score", "TNFA_SIGNALING_VIA_NFKB_Score", "DNA_REPAIR_Score")
scoregroups_test <- paste0(gsub(x = genesets_plot_df$Description, pattern = "HALLMARK_", replacement = ""), "_Score")

# plot top score per tumor, compare high grade vs. low grade ------------------------------------------------------
scoregroup_tmp <- "MTORC1_SIGNALING_Score"
wilcox_pvalues_vec <- NULL
ttest_pvalue_vec <- NULL
# for (scoregroup_tmp in "MTORC1_SIGNALING_Score") {
for (scoregroup_tmp in scoregroups_test) {
  ## make plot data
  scores_tmp_df <- scores_df[, c("cluster_name", scoregroup_tmp)]
  colnames(scores_tmp_df) <- c("cluster_name", "score.bycluster")
  scores_bycase_df <- scores_tmp_df %>%
    mutate(sample_id = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
    mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-")) %>%
    mutate(Case = str_split_fixed(string = Aliquot_WU, pattern = "\\-T|\\-N", n = 2)[,1]) %>%
    group_by(Case) %>%
    summarise(score.bysample = max(score.bycluster))
  plot_data_df <- merge(x = clinical_filtered_df,
                        y = scores_bycase_df,
                        by.x = c("Case"), by.y = c("Case"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Tumor_Stage_Pathological %in% c("Stage I", "Stage II"), "Stage I/II", "Stage III/IV")) %>%
    mutate(y_plot = score.bysample)
  # plot_data_df <- plot_data_df[!duplicated(plot_data_df$Case),]
  
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=.1)
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot), 
                             hide.ns = F, method = "wilcox.test")
  p <- p + ggtitle(label = paste0("Max ", scoregroup_tmp, " by tumor"))
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title = element_blank())
  file2write <- paste0(dir_out, scoregroup_tmp, ".MaxBySample.Wilcox.HighvsLowGrade.png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()
  
  ## wilcox test
  test_result <- wilcox.test(plot_data_df$y_plot[plot_data_df$x_plot == "Stage III/IV"], plot_data_df$y_plot[plot_data_df$x_plot == "Stage I/II"])
  wilcox_pvalues_vec <- c(wilcox_pvalues_vec, test_result$p.value)
  ## t test
  test_result <- t.test(plot_data_df$y_plot[plot_data_df$x_plot == "Stage III/IV"], plot_data_df$y_plot[plot_data_df$x_plot == "Stage I/II"])
  ttest_pvalue_vec <- c(ttest_pvalue_vec, test_result$p.value)
}
result_df <- data.frame(name_score = scoregroups_test, 
                        fdr.wilcox = p.adjust(p = wilcox_pvalues_vec, method = "fdr"), pvalue.wilcox = wilcox_pvalues_vec,
                        fdr.ttest = p.adjust(p = ttest_pvalue_vec, method = "fdr"), pvalue.wilcox = ttest_pvalue_vec)

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "wilcox_ttest_results.tsv")
write.table(x = result_df, file = file2write, quote = F, sep = "\t", row.names = F)

# just test EMT, remove top 3 ---------------------------------------------
for (scoregroup_tmp in "EPITHELIAL_MESENCHYMAL_TRANSITION_Score") {
  ## make plot data
  scores_tmp_df <- scores_df[, c("cluster_name", scoregroup_tmp)]
  colnames(scores_tmp_df) <- c("cluster_name", "score.bycluster")
  scores_bycase_df <- scores_tmp_df %>%
    mutate(sample_id = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
    mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-")) %>%
    mutate(Case = str_split_fixed(string = Aliquot_WU, pattern = "\\-T|\\-N", n = 2)[,1]) %>%
    group_by(Case) %>%
    summarise(score.bysample = max(score.bycluster))
  plot_data_df <- merge(x = clinical_filtered_df,
                        y = scores_bycase_df,
                        by.x = c("Case"), by.y = c("Case"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Tumor_Stage_Pathological %in% c("Stage I", "Stage II"), "Stage I/II", "Stage III/IV")) %>%
    mutate(y_plot = score.bysample) %>%
    filter(!is.na(y_plot)) %>%
    arrange(desc(y_plot))
  cases_remove <- head(plot_data_df$Case, 3)
  plot_data_df <- plot_data_df %>%
    filter(!(Case %in% cases_remove))
  
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=.1)
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot), 
                             hide.ns = F, method = "wilcox.test")
  p <- p + ggtitle(label = paste0("Max ", scoregroup_tmp, " by tumor"))
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title = element_blank())
  file2write <- paste0(dir_out, scoregroup_tmp, ".remove_top3.MaxBySample.Wilcox.HighvsLowGrade.png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()
  
  ## wilcox test
  test_result <- wilcox.test(plot_data_df$y_plot[plot_data_df$x_plot == "Stage III/IV"], plot_data_df$y_plot[plot_data_df$x_plot == "Stage I/II"])
  print(test_result$p.value)
  ## t test
  test_result <- t.test(plot_data_df$y_plot[plot_data_df$x_plot == "Stage III/IV"], plot_data_df$y_plot[plot_data_df$x_plot == "Stage I/II"])
  print(test_result$p.value)
}

# just test EMT, remove 2 outliers ---------------------------------------------
for (scoregroup_tmp in "EPITHELIAL_MESENCHYMAL_TRANSITION_Score") {
  ## make plot data
  scores_tmp_df <- scores_df[, c("cluster_name", scoregroup_tmp)]
  colnames(scores_tmp_df) <- c("cluster_name", "score.bycluster")
  scores_bycase_df <- scores_tmp_df %>%
    mutate(sample_id = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
    mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-")) %>%
    mutate(Case = str_split_fixed(string = Aliquot_WU, pattern = "\\-T|\\-N", n = 2)[,1]) %>%
    group_by(Case) %>%
    summarise(score.bysample = max(score.bycluster))
  plot_data_df <- merge(x = clinical_filtered_df,
                        y = scores_bycase_df,
                        by.x = c("Case"), by.y = c("Case"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Tumor_Stage_Pathological %in% c("Stage I", "Stage II"), "Stage I/II", "Stage III/IV")) %>%
    mutate(y_plot = score.bysample) %>%
    filter(!is.na(y_plot)) %>%
    arrange(desc(y_plot))
  cases_remove <- head(plot_data_df$Case, 2)
  plot_data_df <- plot_data_df %>%
    filter(!(Case %in% cases_remove))
  
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=.1)
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot), 
                             hide.ns = F, method = "wilcox.test")
  p <- p + ggtitle(label = paste0("Max ", scoregroup_tmp, " by tumor"))
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title = element_blank())
  file2write <- paste0(dir_out, scoregroup_tmp, ".remove_outlier2.MaxBySample.Wilcox.HighvsLowGrade.png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()
  
  ## wilcox test
  test_result <- wilcox.test(plot_data_df$y_plot[plot_data_df$x_plot == "Stage III/IV"], plot_data_df$y_plot[plot_data_df$x_plot == "Stage I/II"])
  print(test_result$p.value)
  ## t test
  test_result <- t.test(plot_data_df$y_plot[plot_data_df$x_plot == "Stage III/IV"], plot_data_df$y_plot[plot_data_df$x_plot == "Stage I/II"])
  print(test_result$p.value)
}


# # plot average score per tumor, compare high grade vs. low grade ------------------------------------------------------
# scoregroup_tmp <- "MTORC1_SIGNALING_Score"
# # for (scoregroup_tmp in "MTORC1_SIGNALING_Score") {
# for (scoregroup_tmp in scoregroups_test) {
#   ## make plot data
#   scores_tmp_df <- scores_df[, c("cluster_name", scoregroup_tmp)]
#   colnames(scores_tmp_df) <- c("cluster_name", "score.bycluster")
#   scores_bycase_df <- scores_tmp_df %>%
#     mutate(sample_id = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
#     mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-")) %>%
#     group_by(Aliquot_WU) %>%
#     summarise(score.bysample = mean(score.bycluster))
#   plot_data_df <- merge(x = clinical_filtered_df,
#                         y = scores_bycase_df,
#                         by.x = c("Aliquot.snRNA.WU"), by.y = c("Aliquot_WU"), all.x = T)
#   plot_data_df <- plot_data_df %>%
#     mutate(x_plot = ifelse(Histologic_Grade %in% c("G1", "G2"), "G1/2", "G3/4")) %>%
#     mutate(y_plot = score.bysample) %>%
#     arrange(Case, desc(Histologic_Grade), desc(score.bysample))
#   # plot_data_df <- plot_data_df[!duplicated(plot_data_df$Case),]
#   
#   ## plot
#   p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
#   p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
#   p = p + geom_boxplot(width=.1)
#   p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
#   p = p + stat_compare_means(data = plot_data_df, 
#                              mapping = aes(x = x_plot, y = y_plot), 
#                              hide.ns = F, method = "wilcox.test")
#   p <- p + ggtitle(label = paste0("Mean ", scoregroup_tmp, " by tumor"))
#   p <- p + theme_classic()
#   p <- p + theme(legend.position = "none")
#   p <- p + theme(axis.title = element_blank())
#   file2write <- paste0(dir_out, scoregroup_tmp, ".MeanBySample.Wilcox.HighvsLowGrad.png")
#   png(file2write, width = 600, height = 600, res = 150)
#   print(p)
#   dev.off()
# }
# 
# 
# 
# 
# 
# 
# 
# 
