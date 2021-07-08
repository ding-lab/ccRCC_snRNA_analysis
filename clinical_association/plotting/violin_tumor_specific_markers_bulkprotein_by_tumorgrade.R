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
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_ccRCC_specimen_clinical_data/20210706.v1/ccRCC_Specimen_Clinicl_Data.20210706.v1.tsv")
## input cell type fraction
exp_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Protein/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv")
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# set plotting parameters -------------------------------------------------
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.2, seed = 1)

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

# plot by seperating G1 and G2 --------------------------------------------
# my_comparisons <- list(c("G1", "G3"),c("G1", "G4"), c("G2", "G3"), c("G2", "G4"))
# for (gene_test in genes_process) {
#   ## filter specific protein data
#   exp_test_wide_df <- exp_data_df[exp_df$Index == gene_test,]
#   plot_data_df <- data.frame(CASE_ID = metadata_filtered_df$CASE_ID, Expression = unlist(exp_test_wide_df))
#   plot_data_df$Histologic_Grade <- mapvalues(x = plot_data_df$CASE_ID, from = specimen_clinical_filtered_df$Case, to = as.vector(specimen_clinical_filtered_df$Histologic_Grade))
#   
#   # plot_data_df$easy_id <- barcode_mapping_df$easy_id
#   plot_data_df <- plot_data_df %>%
#     mutate(x_plot = Histologic_Grade) %>%
#     mutate(y_plot = Expression)
#   ## plot
#   p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
#   p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
#   p = p + geom_boxplot(width=.1, alpha = 0.4)
#   p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 1)
#   p = p + stat_compare_means(data = plot_data_df, 
#                              mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), comparisons = my_comparisons,
#                              symnum.args = symnum.args,
#                              hide.ns = F, method = "wilcox.test")
#   p <- p + theme_classic()
#   p <- p + theme(legend.position = "none")
#   p <- p + theme(axis.title.x = element_blank())
#   p <- p + ylab("Bulk protein abundance")
#   file2write <- paste0(dir_out, gene_test, ".png")
#   png(file2write, width = 600, height = 600, res = 150)
#   print(p)
#   dev.off()
# }

# plot by cell group ------------------------------------------------------
my_comparisons <- list(c("G1/2", "G3"),c("G3", "G4"),c("G1/2", "G4"))
for (gene_test in "CP") {
# for (gene_test in genes_process) {
  ## filter specific protein data
  exp_test_wide_df <- exp_data_df[exp_df$Index == gene_test,]
  plot_data_df <- data.frame(CASE_ID = metadata_filtered_df$CASE_ID, Expression = unlist(exp_test_wide_df))
  plot_data_df$Histologic_Grade <- mapvalues(x = plot_data_df$CASE_ID, from = specimen_clinical_filtered_df$Case, to = as.vector(specimen_clinical_filtered_df$Histologic_Grade))
  
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Histologic_Grade %in% c("G1", "G2"), "G1/2", Histologic_Grade)) %>%
    mutate(y_plot = Expression)
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=0.2, alpha = 0.6, outlier.size = 0.1)
  # p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 1)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), comparisons = my_comparisons,
                             symnum.args = symnum.args,
                             hide.ns = F, method = "wilcox.test")
  p <- p + theme_classic(base_size = 12)
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 9))
  p <- p + ylab("Protein level (bulk proteomics)")
  # file2write <- paste0(dir_out, gene_test, ".G12combined.png")
  # png(file2write, width = 600, height = 600, res = 150)
  # print(p)
  # dev.off()
  
  file2write <- paste0(dir_out, gene_test, ".G12combined.pdf")
  pdf(file2write, width = 2, height = 2, useDingbats = F)
  print(p)
  dev.off()
}





