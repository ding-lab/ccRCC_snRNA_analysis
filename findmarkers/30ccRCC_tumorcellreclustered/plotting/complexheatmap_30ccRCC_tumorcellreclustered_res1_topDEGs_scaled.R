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
degs_df <- fread(data.table = F,input = "./Resources/Analysis_Results/findmarkers/30ccRCC_tumorcellreclustered/findallmarker_30ccRCC_tumorcellreclustered_byres/res.1.tumorcellsreclustered.markers.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")

# make matrix data for heatmap body color-----------------------------------------------------------------
## get gene to plot
degs_filtered_df <- degs_df %>%
  filter(p_val_adj < 0.05) %>%
  filter(pct.1 > pct.2) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
degs_filtered_df2 <- degs_df %>%
  filter(gene %in% degs_filtered_df$gene) %>%
  group_by(gene) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
  
genes_plot <- unique(degs_filtered_df2$gene); genes_plot
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
heatmap_colors = circlize::colorRamp2(breaks = c(-1.5, 0, 1.5), colors = c("purple", "black", "yellow"))
heatmap_colors = circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("purple", "black", "yellow"))

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = row_ids, from = degs_filtered_df2$gene, to = as.vector(degs_filtered_df2$cluster))

# make column split -------------------------------------------------------
column_split_vec <- mapvalues(x = column_labels, 
                              from = as.character(c(11, 14, 
                                                    10, 15, 16,
                                                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,12, 13, 17)),
                              to = as.character(c("Cycling", "Cycling", 
                                                  "EMT", "EMT", "EMT",
                                                  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,  12, 13, 17)))

# plot --------------------------------------------------------------------
p <- Heatmap(matrix = plotdata_mat, 
        col = heatmap_colors, 
        cluster_columns = T, show_column_names = T, column_labels = column_labels, column_split = column_split_vec,
        cluster_rows = T, row_split = row_split_vec, cluster_column_slices = T,
        row_names_side = "left", row_dend_side = "right")


# save output -------------------------------------------------------------
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 4
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "heatmap", ".png")
png(file2write, width = 1200, height = 2000, res = 150)
draw(object = p)
dev.off()
