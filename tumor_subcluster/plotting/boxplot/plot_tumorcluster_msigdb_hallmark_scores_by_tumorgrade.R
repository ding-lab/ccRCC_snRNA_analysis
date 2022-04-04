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
# specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_specimen_clinical_data/20200717.v1/snRNA_ccRCC_Specimen_Clinicl_Data.20200717.v1.tsv")
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_ccRCC_specimen_clinical_data/20210706.v1/ccRCC_Specimen_Clinicl_Data.20210706.v1.tsv")
## input cell type fraction
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")

# preprocess -------------------------------------------------
pos <- position_jitter(width = 0.2, seed = 1)
clinical_filtered_df <- specimen_clinical_df %>% 
  filter(Sample_Type == "Tumor") %>%
  filter(Case != "C3L-00359") %>%
  filter(snRNA_available == TRUE) %>%
  select(Aliquot.snRNA.WU, Histologic_Grade, Case) 
scoregroups_test <- colnames(scores_df)
scoregroups_test <- scoregroups_test[grepl(pattern = "Score", x = scoregroups_test)]

# plot top score per tumor, compare high grade vs. low grade ------------------------------------------------------
scoregroup_tmp <- "MTORC1_SIGNALING_Score"

# for (scoregroup_tmp in "MTORC1_SIGNALING_Score") {
for (scoregroup_tmp in scoregroups_test) {
  ## make plot data
  scores_tmp_df <- scores_df[, c("cluster_name", scoregroup_tmp)]
  colnames(scores_tmp_df) <- c("cluster_name", "score.bycluster")
  scores_avg_tmp_df <- scores_tmp_df %>%
    mutate(sample_id = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
    mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-")) %>%
    group_by(Aliquot_WU) %>%
    summarise(score.bysample = max(score.bycluster))
  plot_data_df <- merge(x = clinical_filtered_df,
                        y = scores_avg_tmp_df,
                        by.x = c("Aliquot.snRNA.WU"), by.y = c("Aliquot_WU"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Histologic_Grade %in% c("G1", "G2"), "G1/2", "G3/4")) %>%
    mutate(y_plot = score.bysample) %>%
    arrange(Case, desc(Histologic_Grade), desc(score.bysample))
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
}

# plot average score per tumor, compare high grade vs. low grade ------------------------------------------------------
scoregroup_tmp <- "MTORC1_SIGNALING_Score"
# for (scoregroup_tmp in "MTORC1_SIGNALING_Score") {
for (scoregroup_tmp in scoregroups_test) {
  ## make plot data
  scores_tmp_df <- scores_df[, c("cluster_name", scoregroup_tmp)]
  colnames(scores_tmp_df) <- c("cluster_name", "score.bycluster")
  scores_avg_tmp_df <- scores_tmp_df %>%
    mutate(sample_id = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
    mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-")) %>%
    group_by(Aliquot_WU) %>%
    summarise(score.bysample = mean(score.bycluster))
  plot_data_df <- merge(x = clinical_filtered_df,
                        y = scores_avg_tmp_df,
                        by.x = c("Aliquot.snRNA.WU"), by.y = c("Aliquot_WU"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Histologic_Grade %in% c("G1", "G2"), "G1/2", "G3/4")) %>%
    mutate(y_plot = score.bysample) %>%
    arrange(Case, desc(Histologic_Grade), desc(score.bysample))
  # plot_data_df <- plot_data_df[!duplicated(plot_data_df$Case),]
  
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=.1)
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot), 
                             hide.ns = F, method = "wilcox.test")
  p <- p + ggtitle(label = paste0("Mean ", scoregroup_tmp, " by tumor"))
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title = element_blank())
  file2write <- paste0(dir_out, scoregroup_tmp, ".MeanBySample.Wilcox.HighvsLowGrad.png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()
}








