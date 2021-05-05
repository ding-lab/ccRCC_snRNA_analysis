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
## input mutation info
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/merge_bulk_sn_profiles/20210504.v1/bulk_sn_omics_profile.20210504.v1.tsv")
## input cell type fraction
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/summarize_top_msigdb_geneset_score_by_tumor/20210419.v1/Max.GenesetScorePerSample.tsv")

# set plotting parameters -------------------------------------------------
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.25, seed = 2)
scores_df <- scores_df %>%
  mutate(score_group = variable) %>%
  mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-"))

for (gene_mut in c("SETD2", "KDM5C")) {
  for (scoregroup_tmp in unique(scores_df$score_group)) {
    ## make plot data
    plot_data_df <- mut_df[, c("Aliquot_snRNA_WU", paste0("Mut.", gene_mut))]
    colnames(plot_data_df) <- c("Aliquot_snRNA_WU", "mutation_class")
    plot_data_df <- merge(x = plot_data_df,
                          y = scores_df %>%
                            filter(score_group == scoregroup_tmp) %>%
                            select(Aliquot_WU, value),
                          by.x = c("Aliquot_snRNA_WU"), by.y = c("Aliquot_WU"), all.x = T)
    plot_data_df <- plot_data_df %>%
      filter(!is.na(mutation_class) & !is.na(value)) %>%
      mutate(sample_group = ifelse(mutation_class == "None", "WT", "Mutated")) %>%
      mutate(y_plot = value)
    plot_data_df$x_plot <- factor(x = plot_data_df$sample_group, levels = c("WT", "Mutated"))
    ## plot
    p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
    p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 0.6)
    p = p + geom_boxplot(width=.1, outlier.shape = 23, outlier.fill = "red")
    p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 2)
    p = p + stat_compare_means(data = plot_data_df, 
                               mapping = aes(x = x_plot, y = y_plot), ref.group = "WT",
                               symnum.args = symnum.args,
                               hide.ns = F, method = "wilcox.test")
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





