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
## input
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/avgexp_integrated_data_bycluster_res1_katmai/20220411.v1/30ccRCCtumorcellreclustered.avgexp.RNA.data.byres1clusters.20220411.v1.tsv")
## input the annotation for the hallmark gene sets
degs_df <- fread(data.table = F,input = "./Resources/Analysis_Results/findmarkers/30ccRCC_tumorcellreclustered/findmarker_30ccRCC_tumorcellreclustered_res1_specified_clusters/0_3_7_12_vs_others.v1/res.1.tumorcellsreclustered.markers.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.0_3_7_12_vs_others.v1.tsv")
degs_df <- fread(data.table = F,input = "./Resources/Analysis_Results/findmarkers/30ccRCC_tumorcellreclustered/findmarker_30ccRCC_tumorcellreclustered_res1_specified_clusters/0_3_7_12_vs_others.v2/res.1.tumorcellsreclustered.markers.logfcthreshold.0.minpct.0.mindiffpct.0.0_3_7_12_vs_others.v2.tsv")

# make matrix data for heatmap body color-----------------------------------------------------------------
degs_filtered_df <- degs_df %>%
  mutate(deg_direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
  filter(deg_direction == "up") %>%
  filter(gene_symbol %in% geneset2genes_df$gene[geneset2genes_df$term == "HALLMARK_UV_RESPONSE_DN"]) %>%
  filter(p_val_adj < 0.05)
genes_plot <- degs_filtered_df$gene_symbol; genes_plot
## extract the data for the matrix
plotdata_wide_df <- exp_df %>%
  filter(V1 %in% genes_plot)
plotdata_raw_mat <- as.matrix(plotdata_wide_df[,-1])
## scale across clusters
plotdata_mat <- t(apply(plotdata_raw_mat, 1, scale))
rownames(plotdata_mat) <- plotdata_wide_df$V1
colnames(plotdata_mat) <- colnames(plotdata_raw_mat)

row_ids <- rownames(plotdata_mat)
column_ids <- colnames(plotdata_mat)
column_labels <- gsub(x = column_ids, replacement = "", pattern = "RNA\\.")

# make colors -------------------------------------------------------------
## make colors for the heatmap body
summary(as.vector(plotdata_mat))
heatmap_colors = circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("purple", "black", "yellow"))
heatmap_colors = circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("blue", "white", "red"))

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = row_ids, from = degs_filtered_df$gene_symbol, to = as.vector(degs_filtered_df$deg_direction))

# make column split -------------------------------------------------------
column_split_vec <- mapvalues(x = column_labels, 
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
column_split_vec <- ifelse(column_split_vec == "UV response", "UV response", "other")
# plot --------------------------------------------------------------------
p <- Heatmap(matrix = plotdata_mat, 
        col = heatmap_colors, 
        cluster_columns = T, show_column_names = T, column_labels = column_labels, column_split = column_split_vec,
        cluster_rows = T, cluster_column_slices = F,
        row_split = row_split_vec,
        row_names_side = "left", row_dend_side = "right")


# save output -------------------------------------------------------------
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 5
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "heatmap", ".pdf")
pdf(file2write, width = 6, height = 12, useDingbats = F)
draw(object = p)
dev.off()
file2write <- paste0(dir_out, "heatmap", ".png")
png(file2write, width = 1200, height = 2000, res = 150)
draw(object = p)
dev.off()
