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
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_specimen_clinical_data/20200717.v1/snRNA_ccRCC_Specimen_Clinicl_Data.20200717.v1.tsv")
## input cell type fraction
exp_df <- readRDS(file = "./Resources/Analysis_Results/fetch_data/extract_exp/extract_tumor_specific_marker_sct_data_bybarcode/20210702.v1/tumor_specific_marker_expression_bybarcode.RDS")
barcode_mapping_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/extract_exp/extract_tumor_specific_marker_sct_data_bybarcode/20210702.v1/barcde_mapping.tsv")

# set plotting parameters -------------------------------------------------
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.2, seed = 1)

# preprocess --------------------------------------------------------------
barcode_mapping_df <- barcode_mapping_df %>%
  mutate(aliquot = str_split_fixed(string = aliquot_barcode, pattern = "_", n = 2)[,1])
barcode_mapping_df$Histologic_Grade <- mapvalues(x = barcode_mapping_df$aliquot, from = specimen_clinical_df$Aliquot.snRNA, to = as.vector(specimen_clinical_df$Histologic_Grade))
barcode_mapping_df$easy_id <- mapvalues(x = barcode_mapping_df$aliquot, from = specimen_clinical_df$Aliquot.snRNA, to = as.vector(specimen_clinical_df$Aliquot.snRNA.WU))
barcode_mapping_df <- barcode_mapping_df %>%
  arrange(barcode_merged)
exp_df <- exp_df[barcode_mapping_df$barcode_merged,]

# # plot by seperating G1 and G2 --------------------------------------------
# my_comparisons <- list(c("G1", "G3"),c("G1", "G4"), c("G2", "G3"), c("G2", "G4"))
# for (gene_tmp in colnames(exp_df)) {
#   ## make plot data
#   plot_data_df <- data.frame(barcode_merged = rownames(exp_df), expression = exp_df[,gene_tmp])
#   plot_data_df$Histologic_Grade <- barcode_mapping_df$Histologic_Grade
#   # plot_data_df$easy_id <- barcode_mapping_df$easy_id
#   plot_data_df <- plot_data_df %>%
#     mutate(x_plot = Histologic_Grade) %>%
#     mutate(y_plot = expression)
#   ## plot
#   p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
#   p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
#   p = p + geom_boxplot(width=.1, alpha = 0.4)
#   # p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 1)
#   p = p + stat_compare_means(data = plot_data_df, 
#                              mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), comparisons = my_comparisons,
#                              symnum.args = symnum.args,
#                              hide.ns = F, method = "wilcox.test")
#   p <- p + theme_classic()
#   p <- p + theme(legend.position = "none")
#   p <- p + theme(axis.title.x = element_blank())
#   p <- p + ylab("Tumor-cell Normalized Expression")
#   file2write <- paste0(dir_out, gene_tmp, ".png")
#   png(file2write, width = 600, height = 600, res = 150)
#   print(p)
#   dev.off()
# }

# plot by cell group ------------------------------------------------------
my_comparisons <- list(c("G1/2", "G3"),c("G3", "G4"),c("G1/2", "G4"))
# path_log <- paste0(dir_out, "G12_Combined_FoldChanges.txt")
path_log <- paste0(dir_out, "G12_Combined_T_test.txt")
sink(path_log)
for (gene_tmp in c("UBE2D2", "CP", "NDRG1", "KCTD3", "PCSK6")) {
# for (gene_tmp in colnames(exp_df)) {
  ## make plot data
  plot_data_df <- data.frame(barcode_merged = rownames(exp_df), expression = exp_df[,gene_tmp])
  plot_data_df$Histologic_Grade <- barcode_mapping_df$Histologic_Grade
  plot_data_df$easy_id <- barcode_mapping_df$easy_id
  
  # plot_data_df$easy_id <- barcode_mapping_df$easy_id
  plot_data_df <- plot_data_df %>%
    filter(easy_id != "C3L-00359-T1") %>%
    mutate(x_plot = ifelse(Histologic_Grade %in% c("G1", "G2"), "G1/2", Histologic_Grade)) %>%
    mutate(y_plot = expression)
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=0.2, alpha = 0.4, outlier.size = 0.1)
  # p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 1)
  p = p + stat_compare_means(data = plot_data_df,
                             mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), comparisons = my_comparisons,
                             symnum.args = symnum.args,
                             hide.ns = F, method = "wilcox.test")
  p <- p + ylab("Expression level") + xlab(label = "Tumor grade (number of tumors)")
  p <- p + ggtitle(label = paste0(gene_tmp, " expression "), subtitle = "snRNA-seq")
  p <- p + scale_x_discrete(breaks = c("G1/2", "G3", "G4"), labels = c("G1/2 (12)", "G3 (10)", "G4 (8)"))
  # compute lower and upper whiskers
  p <- p + theme_classic(base_size = 9)
  p <- p + theme(legend.position = "none", title = element_text(size = 9))
  p <- p + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8, color = "black"))
  # file2write <- paste0(dir_out, gene_tmp, ".G12combined.png")
  # png(file2write, width = 600, height = 600, res = 150)
  # print(p)
  # dev.off()
  file2write <- paste0(dir_out, gene_tmp, ".G12combined.pdf")
  pdf(file2write, width = 2, height = 2, useDingbats = F)
  print(p)
  dev.off()
  
  ## calculate fold change
  # exp_med12 <- median(plot_data_df$expression[plot_data_df$Histologic_Grade %in% c("G1", "G2")])
  # fc_3vs12 <- median(plot_data_df$expression[plot_data_df$x_plot == "G3"])/exp_med12
  # fc_4vs12 <- median(plot_data_df$expression[plot_data_df$x_plot == "G4"])/exp_med12
  
  # t.test_result_3vs12 <- t.test(y = plot_data_df$expression[plot_data_df$Histologic_Grade %in% c("G1", "G2")],
  #                            x = plot_data_df$expression[plot_data_df$x_plot == "G3"])
  # t.test_result_4vs12 <- t.test(y = plot_data_df$expression[plot_data_df$Histologic_Grade %in% c("G1", "G2")],
  #                               x = plot_data_df$expression[plot_data_df$x_plot == "G4"])
  # 
  # cat(paste0(gene_tmp, ":\nG3 vs. G1&2: t = ", signif(x = t.test_result_3vs12$statistic, digits = 2), ", p-value = ", signif(x = t.test_result_3vs12$p.value, digits = 2),
  #            "\nG4 vs. G1&2: t = ", signif(x = t.test_result_4vs12$statistic, digits = 2), ", p-value = ", signif(x = t.test_result_4vs12$p.value, digits = 2), "\n\n\n"))
  # cat(paste0(gene_tmp, ":\nG3 vs. G1&2: ", signif(x = fc_3vs12, digits = 2), ";  G4 vs. G1&2: ", signif(x = fc_3vs12, digits = 2), "\n\n\n"))
}
sink()

plot_data_df %>%
  select(easy_id, Histologic_Grade) %>%
  unique() %>%
  select(Histologic_Grade) %>%
  table()



