# Yige Wu @WashU Oct 2020

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
case_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/extract_case_clinical_data/20201125.v1/snRNA_ccRCC_Clinicl_Table.20201125.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input cell type fraction
exp_df <- readRDS(file = "./Resources/Analysis_Results/fetch_data/extract_exp/extract_tumor_specific_marker_sct_data_bybarcode/20210702.v1/tumor_specific_marker_expression_bybarcode.RDS")
barcode_mapping_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/extract_exp/extract_tumor_specific_marker_sct_data_bybarcode/20210702.v1/barcde_mapping.tsv")

# set plotting parameters -------------------------------------------------
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.2, seed = 1)

# preprocess --------------------------------------------------------------
barcode_mapping_df <- barcode_mapping_df %>%
  mutate(aliquot = str_split_fixed(string = aliquot_barcode, pattern = "_", n = 2)[,1])
plot_data_df$Tumor_Stage_Pathological <- mapvalues(x = plot_data_df$Case, from = case_clinical_df$Case, to = as.vector(case_clinical_df$Tumor_Stage_Pathological))
barcode_mapping_df$easy_id <- mapvalues(x = barcode_mapping_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
barcode_mapping_df$Case <- mapvalues(x = barcode_mapping_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
barcode_mapping_df$Tumor_Stage_Pathological <- mapvalues(x = barcode_mapping_df$Case, from = case_clinical_df$Case, to = as.vector(case_clinical_df$Tumor_Stage_Pathological))
barcode_mapping_df <- barcode_mapping_df %>%
  arrange(barcode_merged)
exp_df <- exp_df[barcode_mapping_df$barcode_merged,]

# plot by stage ------------------------------------------------------
my_comparisons <- list(c("Stage I/II", "Stage III"), c("Stage I/II", "Stage IV"), c("Stage III", "Stage IV"))
## make colors for sample group
colors_tumorstage <- brewer.pal(name = "YlOrRd", n = 3)
names(colors_tumorstage) <- c("Stage I/II", "Stage III", "Stage IV")
for (gene_tmp in colnames(exp_df)) {
  ## make plot data
  plot_data_df <- data.frame(barcode_merged = rownames(exp_df), expression = exp_df[,gene_tmp])
  plot_data_df$Tumor_Stage_Pathological <- barcode_mapping_df$Tumor_Stage_Pathological
  # plot_data_df$easy_id <- barcode_mapping_df$easy_id
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Tumor_Stage_Pathological %in% c("Stage I", "Stage II"), "Stage I/II", Tumor_Stage_Pathological)) %>%
    mutate(y_plot = expression)
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=.1, alpha = 0.8)
  p <- p + scale_fill_manual(values = colors_tumorstage)
  # p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 1)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), comparisons = my_comparisons,
                             symnum.args = symnum.args,
                             hide.ns = F, method = "wilcox.test")
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title.x = element_blank())
  p <- p + ylab("Tumor-cell Normalized Expression")
  file2write <- paste0(dir_out, gene_tmp, ".Stage12combined.png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()
}


# plot by seperating stage I & II -----------------------------------------
my_comparisons <- list(c("Stage I", "Stage III"), c("Stage I", "Stage IV"), c("Stage II", "Stage III"), c("Stage II", "Stage IV"))
## make colors for sample group
colors_tumorstage <- brewer.pal(name = "YlOrRd", n = 4)
names(colors_tumorstage) <- c("Stage I", "Stage II", "Stage III", "Stage IV")
for (gene_tmp in colnames(exp_df)) {
  ## make plot data
  plot_data_df <- data.frame(barcode_merged = rownames(exp_df), expression = exp_df[,gene_tmp])
  plot_data_df$Tumor_Stage_Pathological <- barcode_mapping_df$Tumor_Stage_Pathological
  # plot_data_df$easy_id <- barcode_mapping_df$easy_id
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = Tumor_Stage_Pathological) %>%
    mutate(y_plot = expression)
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=.1, alpha = 0.8)
  p <- p + scale_fill_manual(values = colors_tumorstage)
  # p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 1)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), comparisons = my_comparisons,
                             symnum.args = symnum.args,
                             hide.ns = F, method = "wilcox.test")
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title.x = element_blank())
  p <- p + ylab("Tumor-cell Normalized Expression")
  file2write <- paste0(dir_out, gene_tmp, ".png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()
}

