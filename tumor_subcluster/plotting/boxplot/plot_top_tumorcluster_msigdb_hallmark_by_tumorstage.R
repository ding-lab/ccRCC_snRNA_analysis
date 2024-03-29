# Yige Wu @WashU Dec 2020

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
## input clinical data
case_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_case_clinical_data/20201125.v1/snRNA_ccRCC_Clinicl_Table.20201125.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input cell type fraction
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/summarize_top_msigdb_geneset_score_by_tumor/20210929.v1/Max.GenesetScorePerSample.tsv")

# set plotting parameters -------------------------------------------------
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.2, seed = 3)
my_comparisons <- list(c("Stage I/II", "Stage III"), c("Stage I/II", "Stage IV"), c("Stage III", "Stage IV"))
## make colors for sample group
colors_tumorstage <- brewer.pal(name = "YlOrRd", n = 3)
names(colors_tumorstage) <- c("Stage I/II", "Stage III", "Stage IV")

# plot by cell group ------------------------------------------------------
scores_df <- scores_df %>%
  mutate(score_group = variable) %>%
  mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-"))
for (scoregroup_tmp in "EPITHELIAL_MESENCHYMAL_TRANSITION_Score") {
# for (scoregroup_tmp in unique(scores_df$score_gr oup)) {
  ## make plot data
  plot_data_df <- idmetadata_df %>%
    filter(Sample_Type == "Tumor") %>%
    filter(Case != "C3L-00359") %>%
    filter(snRNA_available == TRUE) %>%
    select(Case, Aliquot.snRNA.WU)
  plot_data_df$Tumor_Stage_Pathological <- mapvalues(x = plot_data_df$Case, from = case_clinical_df$Case, to = as.vector(case_clinical_df$Tumor_Stage_Pathological))
  plot_data_df <- merge(x = plot_data_df,
                        y = scores_df %>%
                          filter(score_group == scoregroup_tmp) %>%
                          select(Aliquot_WU, value),
                        by.x = c("Aliquot.snRNA.WU"), by.y = c("Aliquot_WU"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Tumor_Stage_Pathological %in% c("Stage I", "Stage II"), "Stage I/II", Tumor_Stage_Pathological)) %>%
    mutate(y_plot = value) %>%
    arrange(Case, desc(value))
  plot_data_df <- plot_data_df[!duplicated(plot_data_df$Case),]
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p <- p + scale_fill_manual(values = colors_tumorstage)
  p = p + geom_boxplot(width=.1, outlier.shape = 23, outlier.fill = "black")
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 1.5)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), comparisons = my_comparisons,
                             symnum.args = symnum.args,
                             hide.ns = F, method = "wilcox.test")
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  # p <- p + ylab(label = scoregroup_tmp)
  p <- p + ylab(label = "EMT score")
  p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10))
  
  # file2write <- paste0(dir_out, gsub(x = scoregroup_tmp, pattern = "\\/", replacement = "_"), ".png")
  # png(file2write, width = 600, height = 600, res = 150)
  # print(p)
  # dev.off()
  
  file2write <- paste0(dir_out, gsub(x = scoregroup_tmp, pattern = "\\/", replacement = "_"), ".pdf")
  pdf(file2write, width = 3, height = 2.5, useDingbats = F)
  print(p)
  dev.off()
}


# plot --------------------------------------------------------------------
for (scoregroup_tmp in "EPITHELIAL_MESENCHYMAL_TRANSITION_Score") {
  # for (scoregroup_tmp in unique(scores_df$score_gr oup)) {
  ## make plot data
  plot_data_df <- idmetadata_df %>%
    filter(Sample_Type == "Tumor") %>%
    filter(Case != "C3L-00359") %>%
    filter(snRNA_available == TRUE) %>%
    select(Case, Aliquot.snRNA.WU)
  plot_data_df$Tumor_Stage_Pathological <- mapvalues(x = plot_data_df$Case, from = case_clinical_df$Case, to = as.vector(case_clinical_df$Tumor_Stage_Pathological))
  plot_data_df <- merge(x = plot_data_df,
                        y = scores_df %>%
                          filter(score_group == scoregroup_tmp) %>%
                          select(Aliquot_WU, value),
                        by.x = c("Aliquot.snRNA.WU"), by.y = c("Aliquot_WU"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Tumor_Stage_Pathological %in% c("Stage I", "Stage II"), "Stage I/II", "Stage III/IV")) %>%
    mutate(y_plot = value) %>%
    arrange(Case, desc(value))
  plot_data_df <- plot_data_df[!duplicated(plot_data_df$Case),]
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  # p <- p + scale_fill_manual(values = colors_tumorstage)
  p = p + geom_boxplot(width=.1, outlier.shape = 23, outlier.fill = "black")
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 1.5)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), ref.group = "Stage I/II",
                             symnum.args = symnum.args,
                             hide.ns = F, method = "wilcox.test")
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  # p <- p + ylab(label = scoregroup_tmp)
  p <- p + ylab(label = "EMT score")
  p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10), axis.text = element_text(color = "black"))
  
  # file2write <- paste0(dir_out, gsub(x = scoregroup_tmp, pattern = "\\/", replacement = "_"), ".png")
  # png(file2write, width = 600, height = 600, res = 150)
  # print(p)
  # dev.off()
  
  file2write <- paste0(dir_out, gsub(x = scoregroup_tmp, pattern = "\\/", replacement = "_"), ".2groups.", ".pdf")
  pdf(file2write, width = 3, height = 2.5, useDingbats = F)
  print(p)
  dev.off()
}




