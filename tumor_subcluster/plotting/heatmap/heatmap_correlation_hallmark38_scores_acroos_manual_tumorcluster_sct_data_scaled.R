# Yige Wu @WashU Jun 2020

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
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer"
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

# input dependencies ------------------------------------------------------
## input meta data
id_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input pathway scores
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")
## input the pathways to plot
genesets_plot_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/pathway/count_ora_sig_genesets_in_up_degs_across_samples/20220606.v1/Count_gene_set_in_up_tumorcluster_degs.20220606.v1.tsv")
## input the annotation for the hallmark gene sets
hallmark_anno_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Databases/MSigDB/Hallmark_gene_sets_summary.xlsx")
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210805.v1/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

# preprocess --------------------------------------------------------------
genesets_plot_df <- genesets_plot_df %>%
  mutate(scoregroup_name = paste0(gsub(x = Description, pattern = "HALLMARK_", replacement = ""), "_Score")) %>%
  filter(Description != "SPERMATOGENESIS")
genesets_plot <- genesets_plot_df$scoregroup_name
hallmark_anno_df$`Hallmark Name`[hallmark_anno_df$`Hallmark Name` == "UV_RESPONSE_DOWN"] <- "UV_RESPONSE_DN"
## make cluster id
clustername_df <- data.frame(cluster_name = scores_df$cluster_name)
clustername_df <- clustername_df %>%
  mutate(sampleid = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
  mutate(sampleid = gsub(x = sampleid, pattern = "\\.", replacement = "-"))
clustername_df$case <- mapvalues(x = clustername_df$sampleid, from = id_metadata_df$Aliquot.snRNA.WU, to = as.vector(id_metadata_df$Case))
clustername_df <- clustername_df %>%
  filter(case != "C3L-00359") %>%
  filter(!(cluster_name %in% c("C3N.00733.T2_C5", "C3L.01313.T1_C7" , "C3L.01287.T1_C2")))
##
rownames(scores_df) <- scores_df$cluster_name
## prepare CNV data
cnv_df <- cnv_df %>%
  mutate(tumor_subcluster.dataname = gsub(x = tumor_subcluster, pattern = "\\-", replacement = "."))
# format expression data --------------------------------------------------
## get dim names
plot_data_t_mat <- as.matrix(scores_df[,genesets_plot])
plot_data_mat <- t(plot_data_t_mat)
colnames(plot_data_mat) <- scores_df$cluster_name
## make row label
rownames_plot <- rownames(plot_data_mat)
rowlabels_plot <- gsub(x = rownames_plot, pattern = "_Score", replacement = "")
# rowlabels_plot[rowlabels_plot == "EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "EMT"

# make column order -------------------------------------------------------
# col_order_df <- enrich_df %>%
#   arrange(desc(Cell_cycle), desc(Immune), desc(EMT), desc(mTOR))
# plot_data_mat <- plot_data_mat[, col_order_df$cluster_name]
plot_data_mat <- plot_data_mat[, clustername_df$cluster_name]
colnames_plot <- colnames(plot_data_mat)
clusternames_column <- gsub(x = colnames_plot, pattern = "\\.", replacement = "-")
sampleids_column <- str_split_fixed(string = clusternames_column, pattern = "_", n = 2)[,1]

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## colors for gene set scores
summary(as.vector(unlist(plot_data_mat)))
quantile(as.vector(unlist(plot_data_mat)), 0.95)
quantile(as.vector(unlist(plot_data_mat)), 0.05)

colors_correlation <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red")) 

# make row annotation -----------------------------------------------------
freq_de_vec <- mapvalues(x = rownames_plot, from = genesets_plot_df$scoregroup_name, to = as.vector(genesets_plot_df$Freq)); freq_de_vec <- as.numeric(freq_de_vec)
rowanno_obj <- rowAnnotation(Freq_of_DE = anno_barplot(freq_de_vec), annotation_name_side = "top")


# make gene set split -----------------------------------------------------
row_split_vec <- mapvalues(x = rowlabels_plot, from = hallmark_anno_df$`Hallmark Name`, to = as.vector(hallmark_anno_df$`Process Category`))

# plot  ------------------------------------------------------------
p <- Heatmap(matrix = cor(t(plot_data_mat), method = "spearman"),
             col = colors_correlation,
             show_row_names = T, cluster_rows = T,
             cluster_columns = T, 
             # name = "Correlation", 
             show_heatmap_legend = F)

list_lgd = list(
  Legend(col_fun = colors_correlation, 
         title = "Correlation", 
         title_gp = gpar(fontsize = 14),
         labels_gp = gpar(fontsize = 14),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"))
file2write <- paste0(dir_out, "spearman", ".pdf")
pdf(file2write, width = 13, height = 13, useDingbats = F)
draw(object = p,
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

# plot pearson ------------------------------------------------------------
p <- Heatmap(matrix = cor(t(plot_data_mat), method = "pearson"),
             col = colors_correlation,
             show_row_names = T, cluster_rows = T,
             cluster_columns = T, 
             # name = "Correlation", 
             show_heatmap_legend = F)
file2write <- paste0(dir_out, "pearson", ".pdf")
pdf(file2write, width = 13, height = 13, useDingbats = F)
draw(object = p,
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
