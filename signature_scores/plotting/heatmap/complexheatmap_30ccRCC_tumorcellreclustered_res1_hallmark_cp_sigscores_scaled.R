# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

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

rm(sigCorr_df)
# genesets_plot <- sigCorr_df$gene_set[sigCorr_df$FDR < 0.05 & sigCorr_df$C > 0.2]
genesets_plot <- genesets_plot[grepl(pattern = "HALLMARK", x = genesets_plot)]
# genesets_plot <- genesets_plot[grepl(pattern = "HALLMARK|WP_", x = genesets_plot)]
genesets_plot <- genesets_plot[!(genesets_plot %in% c("HALLMARK_UV_RESPONSE"))]
genesets_plot <- c(genesets_plot, "WP_MITOCHONDRIAL_CIV_ASSEMBLY", "HALLMARK_DNA_REPAIR")
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
plotdata_mat <- t(apply(plotdata_raw_mat, 1, scale))
rownames(plotdata_mat) <- plotdata_wide_df$x_plot
colnames(plotdata_mat) <- colnames(plotdata_raw_mat)

row_ids <- rownames(plotdata_mat)
column_ids <- colnames(plotdata_mat)
# row_labels_vec <- gsub(x = row_ids, replacement = "", pattern = "HALLMARK_")
row_labels_vec <- str_split_fixed(string = row_ids, pattern = "_", n = 2)[,2]
# make matrix data for heatmap body significance-----------------------------------------------------------------
## extract the data for the matrix
plotdata_df2 <- results_df %>%
  filter(gene_set %in% genesets_plot) %>%
  mutate(x_plot = gene_set) %>%
  mutate(y_plot = cluster) %>%
  mutate(value = median_sigScore)
plotdata_df3 <- plotdata_df2 %>%
  group_by(gene_set) %>%
  summarise(avg_median_sigScore = mean(value))
plotdata_df2$fdr[plotdata_df2$median_diff < 0] <- 1
plotdata_wide_df2 <- dcast(data = plotdata_df2, formula = x_plot ~ y_plot, value.var = "fdr")
plotdata_mat2 <- as.matrix(plotdata_wide_df2[,-1])
rownames(plotdata_mat2) <- plotdata_wide_df$x_plot
plotdata_mat2 <- plotdata_mat2[row_ids, column_ids]

# make colors -------------------------------------------------------------
## make colors for the heatmap body
summary(as.vector(plotdata_mat))
heatmap_colors = circlize::colorRamp2(breaks = c(-1.5, 0, seq(0.2, 1.8, 0.2)), colors = c("blue", "white", RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")))
## make colors for gene set category
colors_geneset_cat <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")


# make column annotation -----------------------------------------------------
consistency_vec <-  mapvalues(x = row_ids, from = sigCorr_df$gene_set, to = as.vector(sigCorr_df$C)); consistency_vec <- as.numeric(consistency_vec)
geneset_cat_vec <- mapvalues(x = row_labels_vec, from = hallmark_anno_df$`Hallmark Name`, to = as.vector(hallmark_anno_df$`Process Category`))
names(colors_geneset_cat) <- unique(geneset_cat_vec)
anno_obj <- rowAnnotation(Consistency = consistency_vec, Category = anno_simple(x = geneset_cat_vec, col = colors_geneset_cat[geneset_cat_vec]))

# make row split ----------------------------------------------------------
row_split_vec <- geneset_cat_vec

# make column split -------------------------------------------------------
column_split_vec <- mapvalues(x = column_ids, 
                              from = as.character(c(8, 9, 13,
                                                    11, 14, 15,
                                                    0, 3, 7, 12,
                                                    1, 2, 17,
                                                    4, 5, 6, 16,
                                                    10)),
                              to = as.character(c(rep("Immune signaling", 3),
                                                  rep("Proliferating", 3), 
                                                  rep("UV response", 4),
                                                  rep("Metabolic", 3),
                                                  rep("Angiogenic", 4),
                                                  "Unspecified")))

top_n(x = 1:10, n = 5)
# plot --------------------------------------------------------------------
p <- Heatmap(matrix = plotdata_mat, 
        col = heatmap_colors, name = "Scaled signature\nscore",
        cell_fun = function(j, i, x, y, w, h, fill) {
          # if (plotdata_mat[i, j] > 0) {
          #   if(plotdata_mat2[i, j] < 0.001) {
          #     grid.text("***", x, y)
          #   } else if(plotdata_mat2[i, j] < 0.01) {
          #     grid.text("**", x, y)
          #   } else if(plotdata_mat2[i, j] < 0.05) {
          #     grid.text("*", x, y)
          #   }
          # }
          if (plotdata_mat[i,j] >= min(tail(sort(plotdata_mat[,j]), 2))) {
            grid.text("*", x, y)
          }
         },
        cluster_columns = T, show_column_names = T, column_labels = column_split_vec, #column_split = column_split_vec,
        cluster_rows = T, row_split = geneset_cat_vec, cluster_row_slices = T, row_title = NULL,
        #left_annotation = anno_obj, 
        row_names_side = "left", row_labels = row_labels_vec, row_dend_side = "right")

source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "heatmap", ".png")
png(file2write, width = 1500, height = 3000, res = 150)
draw(object = p)
dev.off()
