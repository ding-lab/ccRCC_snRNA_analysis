# Yige Wu @WashU Apr 2020
## plot heatmap the the average expression (not scaled) of HIF pathway members to compare tumor and normal epithelial cells

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

# input denpendencies -----------------------------------------------------
## input average expression by cell type by aliquot
avgexp_bycelltype_byaliquot_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_bycelltypeshorter_byaliquot_on_katmai/20200411.v1/averageexpression_bycelltypeshorter.30_aliquot_integration.20200411.v1.tsv", data.table = F)
## input barcode 2 cell type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200410.v1/30_aliquot_integration.barcode2celltype.20200410.v1.tsv", data.table = F)
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200413.v1/meta_data.20200413.v1.tsv", data.table = F)
## input te bulk genomics/methylation events
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/bulk/other/merge_bulk_sn_profiles/20200512.v1/bulk_sn_omics_profile.20200512.v1.tsv", data.table = F)

# set plot gene list ------------------------------------------------------
gene_plot <- c("VHL")

# make matrix for heatmap body --------------------------------------------
plot_data_df <- avgexp_bycelltype_byaliquot_df %>%
  filter(V1 == gene_plot)
## melt
plot_data_long_df <- melt(plot_data_df)
## add cell type and aliquot info
plot_data_long_df <- plot_data_long_df %>%
  mutate(idaliquot_celltype = gsub(x = variable, pattern = "RNA.", replacement = "")) %>%
  mutate(celltype = str_split_fixed(string = idaliquot_celltype, pattern = "_", n = 2)[,2]) %>%
  mutate(id_aliquot = str_split_fixed(string = idaliquot_celltype, pattern = "_", n = 2)[,1])
plot_data_long_df$Sample_Type <- mapvalues(x = plot_data_long_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Sample_Type))
plot_data_long_df <- plot_data_long_df %>%
  filter(Sample_Type != "Normal")
## dcast
plot_data_wide_df <- dcast(data = plot_data_long_df, formula = celltype ~ id_aliquot, value.var = "value")
## filter columns
rownames(plot_data_wide_df) <- plot_data_wide_df$celltype
data_row_names.keep <- rownames(plot_data_wide_df)
data_row_names.keep <- data_row_names.keep[!(data_row_names.keep %in% c("Unknown"))]
## reformat data frame to matrix
plot_data_mat <- as.matrix(plot_data_wide_df[data_row_names.keep,-1])
plot_data_mat %>% head()
rownames(plot_data_mat) <- plot_data_wide_df[,1]
### get aliquot ids and case ids
ids_aliquot <- colnames(plot_data_mat)
ids_aliquot_wu <- mapvalues(x = ids_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
ids_case <- mapvalues(x = ids_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
ids_celltype <- rownames(plot_data_mat)

# write data matrix -------------------------------------------------------
file2write <- paste0(dir_out, "VHL_input_data_matrix.csv")
write.csv(x = plot_data_mat, file = file2write)

# make column split -------------------------------------------------------
idaliquot_vhl_germline <- bulk_sn_omicsprofile_df$Aliquot.snRNA[!is.na(bulk_sn_omicsprofile_df$Mut.VHL.Germline) & (bulk_sn_omicsprofile_df$Mut.VHL.Germline != "None")]
idaliquot_vhl_germline
idaliquot_3pdel_only <- bulk_sn_omicsprofile_df$Aliquot.snRNA[!is.na(bulk_sn_omicsprofile_df$Mut.VHL) & (bulk_sn_omicsprofile_df$Mut.VHL == "None")]
columnsplit_vec <- ifelse(ids_aliquot %in% idaliquot_vhl_germline, "VHL_Germline",
                          ifelse(ids_aliquot %in% idaliquot_3pdel_only, "3p_Loss", "VHL_Somatic\n+3p_Loss"))
columnsplit_factor <- factor(x = columnsplit_vec, levels = c("VHL_Germline", "VHL_Somatic\n+3p_Loss", "3p_Loss"))

# make heatmap body ------------------------------------------------------------
heatmapbody_color_fun <- colorRamp2(c(quantile(plot_data_mat, 0.025, na.rm=T),
                                      quantile(plot_data_mat, 0.975, na.rm=T),
                                      quantile(plot_data_mat, 1, na.rm=T)),
                                    c("white", 
                                      RColorBrewer::brewer.pal(n = 9, name = "Oranges")[6],
                                      RColorBrewer::brewer.pal(n = 9, name = "Oranges")[9]))

p <- Heatmap(matrix = plot_data_mat,
             col = heatmapbody_color_fun, 
             # cluster_rows = T, row_split = rowsplit_factor, show_row_dend = F,
             # top_annotation = top_col_anno, 
             # show_column_names = F,
             # column_labels = ids_aliquot_wu, 
             column_names_side = "bottom",
             column_split = columnsplit_factor, show_column_dend = F,
             column_title_rot = 90, column_title_side = "top", 
             cluster_rows = F,
             show_column_names = F,
             show_heatmap_legend = F)
p
# make legend -------------------------------------------------------------
## make legend for top annotation
annotation_lgd = list(
  Legend(col_fun = heatmapbody_color_fun, 
         title = "Average expression value\nby cell type (Normalized)", 
         direction = "vertical"))

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "vhl.", "expression.", "tumor_cells.", run_id, ".png")
png(filename = file2write, width = 1500, height = 1200, res = 150)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()

file2write <- paste0(dir_out, "vhl.", "expression.", "tumor_cells.", run_id, ".pdf")
pdf(file = file2write, width = 10, height = 8)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()

