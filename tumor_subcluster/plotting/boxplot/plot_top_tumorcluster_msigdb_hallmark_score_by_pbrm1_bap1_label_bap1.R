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
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210305.v1/meta_data.20210305.v1.tsv")
## input sample group
sample_group_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210304.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210304.v1.tsv")
## input cell type fraction
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/summarize_top_msigdb_geneset_score_by_tumor/20210419.v1/Max.GenesetScorePerSample.tsv")

# set plotting parameters -------------------------------------------------
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.25, seed = 2)
my_comparisons <- list(c("BAP1 mutated", "PBRM1 mutated"), c("BAP1 mutated", "Non-mutants"), c("PBRM1 mutated", "Non-mutants"))
sample_group_df <- sample_group_df %>%
  mutate(sample_group = mutation_category_sim)

# plot by cell group ------------------------------------------------------
scores_df <- scores_df %>%
  mutate(score_group = variable) %>%
  mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-"))
# for (scoregroup_tmp in unique(scores_df$score_group)) {
for (scoregroup_tmp in "MTORC1_SIGNALING_Score") {
    
  ## make plot data
  plot_data_df <- idmetadata_df %>%
    filter(Sample_Type == "Tumor") %>%
    filter(snRNA_available == TRUE) %>%
    select(Case, Aliquot.snRNA.WU)
  plot_data_df$sample_group <- mapvalues(x = plot_data_df$Case, from = sample_group_df$Case, to = as.vector(sample_group_df$sample_group))
  plot_data_df <- merge(x = plot_data_df,
                        y = scores_df %>%
                          filter(score_group == scoregroup_tmp) %>%
                          select(Aliquot_WU, value),
                        by.x = c("Aliquot.snRNA.WU"), by.y = c("Aliquot_WU"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(y_plot = value) %>%
    filter(sample_group != "Both mutated")
  plot_data_df$x_plot <- factor(x = plot_data_df$sample_group, levels = c("Non-mutants", "BAP1 mutated", "PBRM1 mutated", "Both mutated"))
  plot_data_df <- plot_data_df %>%
    mutate(text_label = ifelse(x_plot == "BAP1 mutated", Aliquot.snRNA.WU, NA))
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 0.6)
  p = p + geom_boxplot(width=.1, outlier.shape = 23, outlier.fill = "red")
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), comparisons = my_comparisons,
                             symnum.args = symnum.args,
                             hide.ns = F, method = "wilcox.test")
  p <- p + geom_text_repel(mapping = aes(label = text_label), 
                           max.overlaps = Inf, min.segment.length = unit(0, 'lines'), position = pos)
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 9.5, colour = "black"))
  p <- p + ylab(label = scoregroup_tmp)
  file2write <- paste0(dir_out, gsub(x = scoregroup_tmp, pattern = "\\/", replacement = "_"), ".png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()

}




