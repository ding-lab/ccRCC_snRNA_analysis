# Yige Wu @WashU Aug 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the variable gene list
genes_celltypemarker_df <- fread(input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200831.v1.tsv", data.table = F)
## input the average expression calculated (RNA)
avg.exp.mat <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_bycelltypeshorter_byaliquot_on_katmai/20200828.v1/averageexpression_SCT_bycelltypeshorter_byaliquot.31_aliquot_integration.20200828.v1.tsv", data.table = F)

# specify the genes to show -----------------------------------------------
genes_plot <- genes_celltypemarker_df$Gene[genes_celltypemarker_df$Cell_Type_Group == "Stroma"]

# format the column names to only aliquot id ------------------------------
plot_data_df <- avg.exp.mat %>%
  rename(gene = V1) %>%
  filter(gene %in% genes_plot)
## remove RNA from the column names
data_col_names <- colnames(plot_data_df)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
## rename the data frame
colnames(plot_data_df) <- c("gene", data_col_names.changed)
## remove the subcluster without NA as manual subcluster id
data_col_names.filtered <- data_col_names.changed[grepl(pattern = "_Myofibroblasts", x = data_col_names.changed)]
data_col_names.filtered
plot_data_mat <- as.matrix(plot_data_df[, c(data_col_names.filtered)])
rownames(plot_data_mat) <- plot_data_df$gene
plot_data_mat %>% head()

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
col_fun = colorRamp2(c(quantile(x = plot_data_mat, probs = 0.1, na.rm = T), 
                       quantile(x = plot_data_mat, probs = 0.5, na.rm = T), 
                       quantile(x = plot_data_mat, probs = 0.9, na.rm = T)), 
                     c("white", "yellow", "red"))
## make colors for the discrete ranges
colors_numbercellrange_vec <- RColorBrewer::brewer.pal(n = 6, name = "PuBuGn")
names(colors_numbercellrange_vec) <- sapply(X = seq(from = 1, to = 101, by = 20), FUN = number2rangetext)
## make colors for the cell types
colors_celltype <- stroma_colors

p <- Heatmap(matrix = plot_data_mat, 
             # width = unit(nrow(plot_data_mat), "cm"), height = unit(ncol(plot_data_mat), "cm"),
             # column_labels = ids_aliquot_wu,
             # row_labels = ids_aliquot_wu,
             # right_annotation = row_anno,
             # show_row_names = F, show_column_names = F, 
             # row_split = row_split_factor, cluster_row_slices = F, row_order = row_order_vec,
             # show_row_dend = F, row_title_rot = 0, row_title_side = "right", row_title_gp = gpar(fontsize = 25, fontface = "bold"),
             # row_gap = unit(0, "mm"),
             # top_annotation = col_anno,
             # column_split = column_split_factor, cluster_column_slices = F, column_order = column_order_vec,
             # show_column_dend = F, column_title_side = "top", column_title_rot = 90, column_title_gp = gpar(fontsize = 25, fontface = "bold"),
             # column_gap = unit(0, "mm"),
             # border = "grey50",
             col = col_fun,
             show_row_names = T, show_column_names = T,
             show_heatmap_legend = T)
p
## save heatmap as png
png(filename = paste0(dir_out, "heatmap", ".png"), 
    width = 3000, height = 2400, res = 150)
draw(object = p)
dev.off()


