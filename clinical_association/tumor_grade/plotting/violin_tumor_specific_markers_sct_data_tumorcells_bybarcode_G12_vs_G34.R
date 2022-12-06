# Yige Wu @WashU Oct 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ggpubr"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
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
## input clinical data
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_specimen_clinical_data/20200717.v1/snRNA_ccRCC_Specimen_Clinicl_Data.20200717.v1.tsv")
## input cell type fraction
exp_df <- readRDS(file = "./Resources/Analysis_Results/fetch_data/extract_exp/extract_tumor_specific_marker_sct_data_bybarcode/20221205.v1/tumor_specific_marker_expression_bybarcode.RDS")
barcode_mapping_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/extract_exp/extract_tumor_specific_marker_sct_data_bybarcode/20210702.v1/barcde_mapping.tsv")

# set plotting parameters -------------------------------------------------
my_symnumargs <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.2, seed = 1)
genes_plot <- c("ABCC3", "ABLIM3", "COL23A1", "CP", "EGFR", "ENPP3", "EPHA6", "FTO", "KCTD3", "NDRG1", "PCSK6", "PHKA2", "PLEKHA1", "PLIN2", "SEMA6A", "SHISA9", "SLC6A3", "SNAP25", "TGFA", "UBE2D2" )

# preprocess --------------------------------------------------------------
barcode_mapping_df <- barcode_mapping_df %>%
  mutate(aliquot = str_split_fixed(string = aliquot_barcode, pattern = "_", n = 2)[,1])
barcode_mapping_df$Histologic_Grade <- mapvalues(x = barcode_mapping_df$aliquot, from = specimen_clinical_df$Aliquot.snRNA, to = as.vector(specimen_clinical_df$Histologic_Grade))
barcode_mapping_df$easy_id <- mapvalues(x = barcode_mapping_df$aliquot, from = specimen_clinical_df$Aliquot.snRNA, to = as.vector(specimen_clinical_df$Aliquot.snRNA.WU))
barcode_mapping_df <- barcode_mapping_df %>%
  arrange(barcode_merged)
exp_df <- exp_df[barcode_mapping_df$barcode_merged,]


# plot G1/2 vs. G3/4 ------------------------------------------------------
ttest_t_vec <- NULL
ttest_pvalue_vec <- NULL
for (gene_tmp in genes_plot) {
# for (gene_tmp in colnames(exp_df)) {
  ## make plot data
  plot_data_df <- data.frame(barcode_merged = rownames(exp_df), expression = exp_df[,gene_tmp])
  plot_data_df$Histologic_Grade <- barcode_mapping_df$Histologic_Grade
  plot_data_df$easy_id <- barcode_mapping_df$easy_id
  
  # plot_data_df$easy_id <- barcode_mapping_df$easy_id
  plot_data_df <- plot_data_df %>%
    filter(easy_id != "C3L-00359-T1") %>%
    mutate(x_plot = ifelse(Histologic_Grade %in% c("G1", "G2"), "G1/2", "G3/4")) %>%
    mutate(y_plot = expression)
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  # p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(aes(fill = x_plot), width=0.2, alpha = 0.4, outlier.size = 0.1)
  # p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 1)
  p = p + stat_compare_means(comparisons = list(c("G1/2", "G3/4")), method = "t.test", symnum.args = my_symnumargs)
  p <- p + ylab("Expression level") + xlab(label = "Tumor grade (number of tumors)")
  p <- p + ggtitle(label = paste0(gene_tmp, " expression "), subtitle = "snRNA-seq")
  p <- p + scale_x_discrete(breaks = c("G1/2", "G3/4"), labels = c("G1/2 (12)", "G3/4 (18)"))
  p <- p + ylim(c(0, 1.2*max(plot_data_df$y_plot)))
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
  
  ## do T test
  t.test_result_34vs12 <- t.test(y = plot_data_df$expression[plot_data_df$Histologic_Grade %in% c("G1", "G2")],
                             x = plot_data_df$expression[plot_data_df$Histologic_Grade %in% c("G3", "G4")])

  ttest_pvalue_vec <- c(ttest_pvalue_vec, t.test_result_34vs12$p.value)
  ttest_t_vec <- c(ttest_t_vec, t.test_result_34vs12$statistic)
}

ttest_result_df <- data.frame(gene = genes_plot, p_value = ttest_pvalue_vec, t = ttest_t_vec,
                              fdr = p.adjust(p = ttest_pvalue_vec, method = "fdr"))
file2write <- paste0(dir_out, "t_test.tsv")
write.table(x = ttest_result_df, file = file2write, quote = F, sep = "\t", row.names = F)
