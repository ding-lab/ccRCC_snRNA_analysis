# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/avgeexp_tumorcells_sct_data_by_manualcluster_rm_doublets_on_katmai/20210413.v1/AverageExpression_ByManualTumorSubcluster.20210413.v1.tsv", data.table = F)
## input the genes to plot
pathway2genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/clusterprofiler_tumor_manualsubcluster_up_degs/20210413.v1/ORA_Results.tsv")
## input meta data
id_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210322.v1/meta_data.20210322.v1.tsv")
## input cell number per cluster
cellnumber_percluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/count_cellnumber_per_manual_cluster_rm_doublet/20210413.v1/CellNumberPerTumorManualCluster.20210413.v1.tsv")

-# preprocess --------------------------------------------------------------
## add id for mapping
pathway2genes_df$sample_id <- pathway2genes_df$easy_id
pathway2genes_df$easy_id <- mapvalues(x = pathway2genes_df$sample_id, from = id_metadata_df$Sample, to = as.vector(id_metadata_df$Aliquot.snRNA.WU))
pathway2genes_df <- pathway2genes_df %>%
  mutate(colname_exp = paste0(gsub(x = easy_id,pattern = "\\-", replacement = "."), "_C", (cluster+1)))
## identify clusters with sufficient cell number
cluster_pass_df <- cellnumber_percluster_df %>%
  filter(Freq >= 50)%>%
  mutate(colname_exp = gsub(x = id_cluster_uniq,pattern = "\\-", replacement = "."))

# specify genes to filter -------------------------------------------------
## add name for the marker groups
pathway2genes_filtered_df <- pathway2genes_df %>%
  filter(p.adjust < 0.05) %>%
  filter(grepl(x = Description, pattern = "cell cycle", ignore.case = T))

genes2filter <- sapply(pathway2genes_filtered_df$geneID, function(x) {
  genes_vec <- str_split(string = x, pattern = "\\/")[[1]]
  return(genes_vec)
})
genes2filter <- unique(unlist(genes2filter))

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycluster_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  filter(id_bycluster_byaliquot %in% cluster_pass_df$colname_exp) %>%
  mutate(easyid_column = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(cluster_name = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,2])
## filter out non-tumor and NA tumor cluster
plot_data_long_df <- plot_data_long_df %>%
  filter(!(cluster_name %in% c("", "CNA")))
## make matrix
plot_data_wide_df <- dcast(data = plot_data_long_df, formula = V1 ~ id_bycluster_byaliquot, value.var = "value")
plot_data_raw_mat <- as.matrix(plot_data_wide_df[,-1])
## add row names
rownames(plot_data_raw_mat) <- plot_data_wide_df$V1
## scale by row
plot_data_mat <- t(apply(plot_data_raw_mat, 1, scale))
rownames(plot_data_mat) <- rownames(plot_data_raw_mat)
colnames(plot_data_mat) <- colnames(plot_data_raw_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plot_data_mat))
colors_heatmapbody = colorRamp2(c(-1.5, 
                                  0, 
                                  1.5), 
                                c(color_blue, "white", color_red))
# colors for EMT group
colors_isenriched <- c("black", "grey80")
names(colors_isenriched) <- c("TRUE", "FALSE")

# get row/column ids ----------------------------------------------------------
columnnames_plot <- colnames(plot_data_mat)
ids_aliquot_wu <- str_split_fixed(string = columnnames_plot, pattern = "_", n = 2)[,1]; ids_aliquot_wu
ids_cluster <- str_split_fixed(string = columnnames_plot, pattern = "_", n = 2)[,2]; ids_cluster
rownames_plot <- rownames(plot_data_mat)

# make column annotation --------------------------------------------------
## make cutoff for the EMT-potential-high
cellcycle_enriched_vec <- ifelse(columnnames_plot %in% pathway2genes_filtered_df$colname_exp, "TRUE", "FALSE")
## make column annotation object
colanno_obj = HeatmapAnnotation(CellCycleEnriched = anno_simple(x = cellcycle_enriched_vec, col = colors_isenriched[cellcycle_enriched_vec], height = unit(1, "cm")),
                                annotation_name_gp = gpar(fontsize = 15, fontface = "bold"), annotation_name_side = "left")

# make column order --------------------------------------------------
# column_order_vec <- order(scores_mesenchymal, decreasing = T)
# emt_group_factor <- factor(x = emt_group_vec, levels = c("Low", "High"))
# column_order_vec <- order(emt_group_factor, scores_mesenchymal, decreasing = T)

#column_order_vec

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na, border = "black",
                             ## row
                             show_row_names = T, row_names_gp = gpar(fontsize = 10, fontface = "italic"), row_names_side = "left",
                             # row_split = row_split_factor,
                             row_title_rot = 0, row_title_gp = gpar(fontsize = 20),
                             # row_labels = factor_cellgroup,
                             show_row_dend = F, 
                             ## column
                             show_column_dend = F, cluster_columns = T, 
                             # column_order = column_order_vec,
                             top_annotation = colanno_obj, show_column_names = T, column_title = NA,
                             show_heatmap_legend = F)
p


# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_heatmapbody, 
         title = "snRNA expression (scaled by row)", 
         title_gp = gpar(fontsize = 20),
         labels_gp = gpar(fontsize = 20),
         legend_width = unit(6, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cell_Cycle_Genes_by_tumorcluster", ".png")
png(file2write, width = 2500, height = 1500, res = 150)
draw(object = p, 
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()
# file2write <- paste0(dir_out, "Cell_Cycle_Genes_by_tumorcluster", ".pdf")
# pdf(file2write, width = 20, height = 6.5, useDingbats = F)
# draw(object = p, 
#      annotation_legend_side = "top", annotation_legend_list = list_lgd)
# dev.off()

