# Yige Wu @WashU May 2022

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
## input median(?) signature scores per cluster
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/compare_signature_scores/wilcox_eachcluster_vs_rest_bygeneset_res1_30ccRCC_tumorcellsreclustered/20220411.v2/wilcox.eachcluster_vs_rest.res1.30ccRCCtumorcellreclustered,20220411.v2.tsv")
## input the auto-correlation results
sigCorr_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/run_vision/getSignatureAutocorrelation_30ccRCC_tumorcellreclustered/20220411.v1/ccRCC.30ccRCC.TumorCellsReclustered.Vision.SignatureAutocorrelation.20220411.v1.tsv")
## input the annotation for the hallmark gene sets
hallmark_anno_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Databases/MSigDB/Hallmark_gene_sets_summary.xlsx")

# make matrix data for heatmap body color-----------------------------------------------------------------
## get gene sets to plot
genesets_plot <- sigCorr_df$gene_set[sigCorr_df$FDR < 0.05 & sigCorr_df$C > 0.1]
# genesets_plot <- sigCorr_df$gene_set[sigCorr_df$FDR < 0.05 & sigCorr_df$C > 0.2]
genesets_plot <- genesets_plot[grepl(pattern = "HALLMARK", x = genesets_plot)]
genesets_plot <- genesets_plot[!(genesets_plot %in% c("HALLMARK_UV_RESPONSE"))]
rm(sigCorr_df)
## extract the data for the matrix
plotdata_df <- results_df %>%
  filter(gene_set %in% genesets_plot) %>%
  mutate(x_plot = gene_set) %>%
  mutate(y_plot = cluster) %>%
  mutate(value = median_sigScore)
rm(results_df)
plotdata_wide_df <- dcast(data = plotdata_df, formula = x_plot ~ y_plot, value.var = "value")
plotdata_raw_mat <- as.matrix(plotdata_wide_df[,-1])
## scale across clusters
plotdata_mat <- apply(plotdata_raw_mat, 1, scale)
rownames(plotdata_mat) <- colnames(plotdata_raw_mat)
colnames(plotdata_mat) <- plotdata_wide_df$x_plot
## reorder clusters
clusters_ordered <- rev(c("10", "15", "6", "13", "9", "8", "5", "14", "16", "11", "4", "12", "7", "0", "3", "17", "2", "1"))
plotdata_mat <- plotdata_mat[clusters_ordered,]

# make row labels ---------------------------------------------------------
row_labels_vec <- c("MC1:ER_signaling", "MC2:OXPHOS_high", "MC17:OXOHOS_high",
                "MC3:DDR_dysregulated", "MC0:DDR_dysregulated", "MC7:DDR_dysregulated", "MC12:DDR_dysregulated",
                "MC4:ER_signaling", "MC11:Cycling", "MC16:Multi-activated", "MC14:Cycling", 
                "MC5:IL2-STAT5_signaling", "MC8:IFN_signaling", "MC9:Inflammatory", "MC13:Inflammatory+KRAS signaling",
                "MC6:Hypoxic", "MC15:Multi-activated", "MC10:Multi-activated")
row_labels_vec <- c("MC2:ER_signaling", "MC3:OXPHOS_high", "MC18:OXOHOS_high",
                    "MC4:DDR_dysregulated", "MC1:DDR_dysregulated", "MC8:DDR_dysregulated", "MC13:DDR_dysregulated",
                    "MC5:ER_signaling", "MC12:Cycling", "MC17:Multi-activated", "MC15:Cycling", 
                    "MC6:IL2-STAT5_signaling", "MC9:IFN_signaling", "MC10:Inflammatory", "MC14:Inflammatory+KRAS signaling",
                    "MC7:Hypoxic", "MC16:Multi-activated", "MC11:Multi-activated")

row_ids <- rownames(plotdata_mat)
column_ids <- colnames(plotdata_mat)
col_labels_vec <- gsub(x = column_ids, replacement = "", pattern = "HALLMARK_")
col_labels_vec <- tolower(col_labels_vec)
col_labels_vec[col_labels_vec == tolower("EPITHELIAL_MESENCHYMAL_TRANSITION")] <- "EMT"
col_labels_vec <- gsub(x = col_labels_vec, pattern = "uv_", replacement = "UV_")
col_labels_vec <- gsub(x = col_labels_vec, pattern = "e2f_", replacement = "E2F_")
col_labels_vec <- gsub(x = col_labels_vec, pattern = "g2m_", replacement = "G2M_")
col_labels_vec <- gsub(x = col_labels_vec, pattern = "myc_", replacement = "MYC_")
col_labels_vec <- gsub(x = col_labels_vec, pattern = "il6_jak_stat3", replacement = "IL6-JAK-STAT3")
col_labels_vec <- gsub(x = col_labels_vec, pattern = "il2_stat5", replacement = "IL2-STAT5")
col_labels_vec <- gsub(x = col_labels_vec, pattern = "pi3k_akt_mtor", replacement = "PI3K-AKT-mTOR")
col_labels_vec <- gsub(x = col_labels_vec, pattern = "kras", replacement = "KRAS")
col_labels_vec[col_labels_vec == "tnfa_signaling_via_nfkb"] <- expression("TNF-alpha signaling via NFkB")
col_labels_vec <- gsub(x = col_labels_vec, pattern = "mtorc1", replacement = "mTORC1")
col_labels_vec <- gsub(x = col_labels_vec, pattern = "tgf_beta", replacement = "TGF-beta")
col_labels_vec <- gsub(x = col_labels_vec, pattern = "_", replacement = " ")

# make colors -------------------------------------------------------------
## make colors for the heatmap body
summary(as.vector(plotdata_mat))
# colors_sig_zscore = circlize::colorRamp2(breaks = c(-1.5, 0, seq(0.2, 1.8, 0.2)), 
#                                       colors = c("blue", "white", RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")))
colors_sig_zscore = circlize::colorRamp2(breaks = c(-2, 0, 2), 
                                      colors = c("purple", "black", "yellow"))
## 
colors_correlation <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red")) 

# make row annotation -----------------------------------------------------
# consistency_vec <-  mapvalues(x = row_ids, from = sigCorr_df$gene_set, to = as.vector(sigCorr_df$C)); consistency_vec <- as.numeric(consistency_vec)
geneset_cat_vec <- mapvalues(x =  gsub(x = column_ids, replacement = "", pattern = "HALLMARK_"), from = hallmark_anno_df$`Hallmark Name`, to = as.vector(hallmark_anno_df$`Process Category`))
names(colors_geneset_cat) <- unique(geneset_cat_vec)
# row_anno_obj <- rowAnnotation(Category = anno_simple(x = geneset_cat_vec, col = colors_geneset_cat[geneset_cat_vec]))

# make column split ----------------------------------------------------------
col_split_vec <- geneset_cat_vec
col_split_vec <- factor(x = col_split_vec, levels = c("metabolic", "DNA damage", "proliferation", "immune", "development", "pathway", "signaling", "cellular component"))

# plot --------------------------------------------------------------------
fontsize_plot <- 14
p <- Heatmap(matrix = cor(t(plotdata_mat), method = "spearman"),
                 col = colors_correlation,
             show_row_names = T, cluster_rows = T, show_row_dend = F,
             cluster_columns = T, column_labels = paste0("MC", as.numeric(clusters_ordered)+1), column_names_gp = gpar(fontsize = fontsize_plot), column_dend_side = "bottom",
             # name = "Correlation", 
             show_heatmap_legend = F)
p <- p + Heatmap(matrix = plotdata_mat, 
             col = colors_sig_zscore,
             cell_fun = function(j, i, x, y, w, h, fill) {
               if (plotdata_mat[i,j] >= quantile(plotdata_mat[,j], 0.75)+1.5*IQR(plotdata_mat[,j])) {
                 grid.text("*", x, y)
               }
               if (plotdata_mat[i,j] >= max(plotdata_mat[i,])) {
               # if (plotdata_mat[i,j] >= min(tail(sort(plotdata_mat[i,]), 2))) {
                 grid.rect(x = x, y = y, width = w, height = h, 
                           gp = gpar(col = "red", fill = NA))
               }
             },
             cluster_rows = F,
             show_row_names = T,
             # top_annotation = col_anno_obj,
             cluster_columns = T, show_column_dend = F,
             column_split = col_split_vec, column_title_gp = gpar(fontsize = fontsize_plot), column_title_rot = 90,
             cluster_column_slices = F, 
             column_labels = col_labels_vec, column_names_gp = gpar(fontsize = fontsize_plot),
             # # left_annotation = row_anno_obj, 
             row_names_side = "right", row_labels = row_labels_vec, row_names_gp = gpar(fontsize = fontsize_plot),
             show_heatmap_legend = F)

list_lgd = list(
  Legend(col_fun = colors_correlation, 
         title = "Correlation",
         title_gp = gpar(fontsize = fontsize_plot),
         labels_gp = gpar(fontsize = fontsize_plot),
         # legend_width = unit(4, "cm"),
         # legend_height = unit(4, "cm"),
         direction = "vertical"),
  Legend(col_fun = colors_sig_zscore, 
         title = "Signature\nz-score",
         title_gp = gpar(fontsize = fontsize_plot),
         labels_gp = gpar(fontsize = fontsize_plot),
         # legend_width = unit(4, "cm"),
         # legend_height = unit(4, "cm"),
         direction = "vertical"))


p

source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "heatmap", ".png")
png(file2write, width = 2000, height = 1000, res = 150)
draw(object = p,
     annotation_legend_side = "left", annotation_legend_list = list_lgd)
dev.off()

file2write <- paste0(dir_out, "heatmap", ".pdf")
pdf(file2write, width = 15, height = 7, useDingbats = F)
draw(object = p,
     annotation_legend_side = "left", annotation_legend_list = list_lgd)
dev.off()

