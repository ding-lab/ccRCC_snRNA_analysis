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
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_ccRCC_specimen_clinical_data/20210706.v1/ccRCC_Specimen_Clinicl_Data.20210706.v1.tsv")
## input cell type fraction
exp_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Protein/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv")
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# set plotting parameters -------------------------------------------------
my_symnumargs <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.2, seed = 1)
genes_plot <- c("ABCC3", "ABLIM3", "COL23A1", "CP", "EGFR", "ENPP3", "EPHA6", "FTO", "KCTD3", "NDRG1", "PCSK6", "PHKA2", "PLEKHA1", "PLIN2", "SEMA6A", "SHISA9", "SLC6A3", "SNAP25", "TGFA", "UBE2D2" )

# preprocess --------------------------------------------------------------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma")
## center data
exp_data_df <- exp_df[, metadata_filtered_df$Specimen.Label.tumor]
exp_data_df <- (exp_data_df - exp_df$ReferenceIntensity)
## rename columns
colnames(exp_data_df) <- metadata_filtered_df$CASE_ID
## filter clinical data
specimen_clinical_filtered_df <- specimen_clinical_df %>%
  filter(discovery_study == "Yes") %>%
  filter(Tissue_Type == "tumor")

# specify gene to test ----------------------------------------------------
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210702.v1/ccRCC_markers.Surface.20210702.v1.tsv")
genes_process <- genes_process_df$Gene
genes_process <- genes_process[genes_process %in% exp_df$Index]


# plot G1/2 vs. G3/4 ------------------------------------------------------
ttest_t_vec <- NULL
ttest_pvalue_vec <- NULL
fontsize <- 10

# for (gene_tmp in genes_plot) {
for (gene_tmp in c("CP", "PCSK6")) {
  ## filter specific protein data
  exp_test_wide_df <- exp_data_df[exp_df$Index == gene_tmp,]
  if (nrow(exp_test_wide_df) == 0) {
    ttest_pvalue_vec <- c(ttest_pvalue_vec, NA)
    ttest_t_vec <- c(ttest_t_vec, NA)
    next()
  }
  plot_data_df <- data.frame(CASE_ID = metadata_filtered_df$CASE_ID, Expression = unlist(exp_test_wide_df))
  plot_data_df$Histologic_Grade <- mapvalues(x = plot_data_df$CASE_ID, from = specimen_clinical_filtered_df$Case, to = as.vector(specimen_clinical_filtered_df$Histologic_Grade))
  
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Histologic_Grade %in% c("G1", "G2"), "G1/2", "G3/4")) %>%
    mutate(y_plot = Expression) %>%
    filter(!is.na(y_plot))
  number_tumors_df <- plot_data_df %>%
    select(CASE_ID, x_plot) %>%
    unique() %>%
    select(x_plot) %>%
    table() %>%
    as.data.frame() %>%
    rename(x_plot = '.')
  number_G12 <- number_tumors_df$Freq[number_tumors_df$x_plot == "G1/2"]; x_G12 <- paste0("G1/2 (", number_G12, ")")
  number_G34 <- number_tumors_df$Freq[number_tumors_df$x_plot == "G3/4"]; x_G34 <- paste0("G3/4 (", number_G34, ")")
  plot_data_df <- plot_data_df %>%
    select(x_plot, y_plot)
  
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p <- p + stat_boxplot(geom = "errorbar", width = 0.2)
  p <- p + geom_boxplot(width=0.2, alpha = 0.6, outlier.size = 0.1)
  p = p + stat_compare_means(comparisons = list(c("G1/2", "G3/4")), method = "t.test", symnum.args = my_symnumargs)
  p <- p + ylab("Expression level") + xlab(label = "Tumor grade (No. tumors)")
  p <- p + ggtitle(label = paste0(gene_tmp, " protein "), subtitle = "bulk proteomics")
  p <- p + scale_x_discrete(breaks = c("G1/2", "G3/4"), 
                            labels = c(x_G12, x_G34))
  p <- p + ylim(c(min(plot_data_df$y_plot)-IQR(plot_data_df$y_plot)*0.2, max(plot_data_df$y_plot)+IQR(plot_data_df$y_plot)*0.7))
  p <- p + theme_classic(base_size = fontsize)
  p <- p + theme(legend.position = "none", title = element_text(size = fontsize))
  p <- p + theme(axis.title = element_text(size = fontsize), axis.text = element_text(size = fontsize, color = "black"))
  
  file2write <- paste0(dir_out, gene_tmp, ".pdf")
  pdf(file2write, width = 2, height = 2, useDingbats = F)
  print(p)
  dev.off()
  
  ## write source data
  write.table(x = plot_data_df, file = paste0("~/Desktop/SF2d.", gene_tmp, ".SourceData.tsv"), quote = F, sep = "\t", row.names = F)
  
  # ## do T test
  # t.test_result_34vs12 <- t.test(y = plot_data_df$y_plot[plot_data_df$Histologic_Grade %in% c("G1", "G2")],
  #                                x = plot_data_df$y_plot[plot_data_df$Histologic_Grade %in% c("G3", "G4")])
  # 
  # ttest_pvalue_vec <- c(ttest_pvalue_vec, t.test_result_34vs12$p.value)
  # ttest_t_vec <- c(ttest_t_vec, t.test_result_34vs12$statistic)
}

# ttest_result_df <- data.frame(gene = genes_plot, p_value = ttest_pvalue_vec, t = ttest_t_vec,
#                               fdr = p.adjust(p = ttest_pvalue_vec, method = "fdr"))
# file2write <- paste0(dir_out, "t_test.tsv")
# write.table(x = ttest_result_df, file = file2write, quote = F, sep = "\t", row.names = F)
# 







