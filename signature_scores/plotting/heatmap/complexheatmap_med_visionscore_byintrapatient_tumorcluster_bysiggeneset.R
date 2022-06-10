
# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "circlize",
  "ComplexHeatmap"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input the auto-correlation results
sigCorr_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/run_vision/getSignatureAutocorrelation_30ccRCC_tumorcellreclustered/20220411.v1/ccRCC.30ccRCC.TumorCellsReclustered.Vision.SignatureAutocorrelation.20220411.v1.tsv")
## input the annotation for the hallmark gene sets
hallmark_anno_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Databases/MSigDB/Hallmark_gene_sets_summary.xlsx")
## input median(?) signature scores per cluster
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/process_signature_scores/getmedian_visionscore_bysiggeneset_byintrapatientcluster/20220606.v1/median_vision_score.sig_genesets.byintrapatient_tumor_cluster.20220606.v1.tsv")
## input clusters to plot scores
ics_plot_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")
## input barcode count by meta-cluster
barcode_bymc_byic_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/other/count_intrapatient_tumorcluster_cells_by_meta_cluster_res1/20220606.v1/median_vision_score.sig_genesets.byintrapatient_tumor_cluster.20220606.v1.tsv")

# pre-process -------------------------------------------------------------
## get gene sets to plot
genesets_plot <- sigCorr_df$gene_set[sigCorr_df$FDR < 0.05 & sigCorr_df$C > 0.1]
# genesets_plot <- sigCorr_df$gene_set[sigCorr_df$FDR < 0.05 & sigCorr_df$C > 0.2]
genesets_plot <- genesets_plot[grepl(pattern = "HALLMARK", x = genesets_plot)]
genesets_plot <- genesets_plot[!(genesets_plot %in% c("HALLMARK_UV_RESPONSE"))]
rm(sigCorr_df)
##
ic_plot <- ics_plot_df$cluster_name[!grepl(pattern = "359", x = ics_plot_df$cluster_name)]
ic_plot <- gsub(x = ic_plot, pattern = "\\.", replacement = "-")

# make plot matrix for signature score --------------------------------------------------------
## extract the data for the matrix
plotdata_raw_mat <- t(as.matrix(results_df[,genesets_plot]))
## scale across clusters
plotdata_mat <- t(apply(plotdata_raw_mat, 1, scale))
rownames(plotdata_mat) <- rownames(plotdata_raw_mat)
colnames(plotdata_mat) <- results_df$intrapatient_cluster_name
plotdata_mat <- plotdata_mat[,ic_plot]
row_ids_vec <- rownames(plotdata_mat)

ic_plot[!(ic_plot %in% results_df$intrapatient_cluster_name)]

# make plot matrix for barcode % --------------------------------------------------------
barcode_bymc_byic_df <- barcode_bymc_byic_df %>%
  mutate(perc_bymc_byic = (cells_bymc_byic/cells_byic)*100)
plotdata_df2 <- dcast(data = barcode_bymc_byic_df, formula = meta_cluster_name ~ intrapatient_cluster_name, value.var = "perc_bymc_byic")
plotdata_df2[is.na(plotdata_df2)] <- 0
plotdata_mat2 <- as.matrix(plotdata_df2[,ic_plot])
rownames(plotdata_mat2) <- plotdata_df2$meta_cluster_name

# make colors -------------------------------------------------------------
## make colors for the heatmap body
summary(as.vector(plotdata_mat))
colors_sig_zscore = circlize::colorRamp2(breaks = c(-2, 0, 2), 
                                         colors = c("purple", "black", "yellow"))


# make split and labels for gene sets ------------------------------------------------
row_labels_vec <- gsub(x = row_ids_vec, replacement = "", pattern = "HALLMARK_")
geneset_cat_vec <- mapvalues(x = row_labels_vec, from = hallmark_anno_df$`Hallmark Name`, to = as.vector(hallmark_anno_df$`Process Category`))
row_split_vec <- geneset_cat_vec
row_split_vec <- factor(x = row_split_vec, levels = c("metabolic", "DNA damage", "proliferation", "immune", "development", "pathway", "signaling", "cellular component"))


# plot --------------------------------------------------------------------
p <- Heatmap(matrix = plotdata_mat, 
             col = colors_sig_zscore,
             # cell_fun = function(j, i, x, y, w, h, fill) {
             #   if (plotdata_mat[i,j] >= quantile(plotdata_mat[,j], 0.75)+1.5*IQR(plotdata_mat[,j])) {
             #     grid.text("*", x, y)
             #   }
             #   if (plotdata_mat[i,j] >= max(plotdata_mat[i,])) {
             #     # if (plotdata_mat[i,j] >= min(tail(sort(plotdata_mat[i,]), 2))) {
             #     grid.rect(x = x, y = y, width = w, height = h, 
             #               gp = gpar(col = "red", fill = NA))
             #   }
             # },
             cluster_rows = F,
             show_row_names = T,
             # top_annotation = col_anno_obj,
             cluster_columns = T,
             cluster_column_slices = F, 
             # # left_annotation = row_anno_obj, 
             row_names_side = "right", row_labels = row_labels_vec, row_split = row_split_vec,
             show_heatmap_legend = F)
p <- Heatmap(matrix = cor(t(plotdata_mat), method = "spearman"),
             col = colors_correlation,
             show_row_names = T, cluster_rows = T,
             cluster_columns = T, column_labels = paste0("MC", clusters_ordered),
             # name = "Correlation", 
             show_heatmap_legend = F)
p




