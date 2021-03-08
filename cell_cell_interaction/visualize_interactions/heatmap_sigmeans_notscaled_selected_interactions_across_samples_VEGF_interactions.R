# Yige Wu @WashU Sep 2020

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
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/rank_across_celltypes_by_pair_by_avgsigmean/20201012.v1/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
## input filtered summary data
interactions_filtered_df <- fread(data.table = F, input = "./Data_Freezes/V1/snRNA/Cell_Cell_Interactions/Values_by_sample_mean_expr_of_curated_interactions_top3ct_stroma_immune_20200818_seen_5+_grouped_by_Inter.group_strongest_in_given_pair.txt")
## input cell-cell interaction values by case
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/other/scale_pair_cell_types_vs_sample_sig_mean/20200924.v1/cellphonedb.pair_cell.types_vs_sample.sig_mean.tsv")
# cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/filter_interactions/filter_cellphonedb_out/20201012.v1/cell.phone.res.total.run20200818.filtered.txt")
## input sample meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input clinical info
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/extract_specimen_clinical_data/20200717.v1/snRNA_ccRCC_Specimen_Clinicl_Data.20200717.v1.tsv")

# specify pairs to filter -------------------------------------------------
## get genes to filter
summary_filtered_df <- summary_df %>%
  mutate(interacting_pair_directed = paste0(gene.source, "_", gene.target)) %>%
  filter(interacting_pair_directed %in% c("VEGFA_FLT1", "VEGFB_FLT1", 'VEGFA_KDR')) %>%
  filter(pair_cell.types %in% interactions_filtered_df$pair_cell.types[interactions_filtered_df$Inter.group != "Immune:Stroma"]) %>%
  # filter(pair_cell.types %in% interactions_filtered_df$pair_cell.types | celltypes.source2target %in% c("Normal epithelial cells->Endothelial cells")) %>%
  arrange(desc(avg_sig_mean))
interacting_pairs_process <- summary_filtered_df$interacting_pair
interacting_pairs_process
length(interacting_pairs_process)
paired_celltypes_process <- summary_filtered_df$paired_celltypes
genes_process <- unique(c(summary_filtered_df$gene.source, summary_filtered_df$gene.target))

# specify the aliquot ids (ordered) ---------------------------------------

# make data matrix for the heatmap body -------------------------------------------
plot_data_df <- cellphone_df %>%
  filter(pair_cell.types %in% summary_filtered_df$pair_cell.types)
plot_data_mat <- as.matrix(plot_data_df[, -1])
rownames(plot_data_mat) <- plot_data_df[,1]
## remove the sample without endothelial cells
plot_data_mat <- plot_data_mat[,!(colnames(plot_data_mat) %in% c("C3L-00359-T1", "C3L-00088-N"))]
# plot_data_mat[1:3, 1:10]
# summary(as.vector(plot_data_mat))
##
plot_data_mat[is.na(plot_data_mat)] <- -3
ids_aliquot <- colnames(plot_data_mat)
interaction_celltypes <- rownames(plot_data_mat)

# preprocess the protein data ---------------------------------------------
## get the bulk aliquot ids
ids_aliquot_bulk <- mapvalues(x = ids_aliquot, from = idmetadata_df$Aliquot.snRNA.WU, to = as.vector(idmetadata_df$Aliquot.bulk))
## make the matrix to plot the heatmap
protein_df2plot <- protein_df %>%
  filter(Index %in% genes_process) %>%
  select("Index", "ReferenceIntensity", ids_aliquot_bulk[!is.na(ids_aliquot_bulk)])

protein_mat1 <- protein_df2plot %>%
  select(-Index) %>%
  select(-ReferenceIntensity)
rownames(protein_mat1) <- protein_df2plot$Index
protein_mat2 <- protein_mat1
protein_mat2 <- as.matrix(protein_mat1) - as.vector(protein_df2plot$ReferenceIntensity)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "white"
## make color function for heatmap body color
RColorBrewer::display.brewer.all()
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
color_orange <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[5]
color_yellow <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")[1]
summary(as.vector(plot_data_mat[plot_data_mat != -3]))
colors_heatmapbody <- colorRamp2(c(-3, 0, 10, 20, 30, 40), 
                                 c(color_na, brewer.pal(n = 5, name = "YlOrRd")))
colors_heatmapbody_legand <- colorRamp2(c(0, 10, 20, 30, 40), 
                                        c(brewer.pal(n = 5, name = "YlOrRd")))
## cell type colors
unique(c(plot_data_df$Cell_type.source, plot_data_df$Cell_type.target))
colors_celltype <- rep(x = colors_cellgroup14[c("Tumor cells", "Endothelial cells", "Myofibroblasts", "Fibroblasts", "Macrophages", "NK cells")], c(1, 1, 1, 1, 3, 2))
names(colors_celltype) <- c("Tumor cells", "Endothelial cells", "Myofibroblasts", "Fibroblasts", 
                            "Macrophages", "Macrophages proliferating", "TRM", 
                            "NK cells strong", "NK cells weak")
## make colors for sample type
colors_sampletype <- c("Tumor" = color_red, "Normal" = RColorBrewer::brewer.pal(n = 3, name = "Set1")[3])

# make row annotation -----------------------------------------------------
avg_sig_mean_vec <- mapvalues(x = interaction_celltypes, from = summary_df$pair_cell.types, to = as.vector(summary_df$avg_sig_mean))
avg_sig_mean_vec <- as.numeric(avg_sig_mean_vec)
avg_sig_mean_vec
summary(avg_sig_mean_vec)
## make text
gene_source_vec <- mapvalues(x = interaction_celltypes, from = summary_filtered_df$pair_cell.types, to = as.vector(summary_filtered_df$gene.source))
gene_target_vec <- mapvalues(x = interaction_celltypes, from = summary_filtered_df$pair_cell.types, to = as.vector(summary_filtered_df$gene.target))
gene_text_vec <- paste0(gene_source_vec, "->", gene_target_vec)
## make cell types
celltype_source_vec <- mapvalues(x = interaction_celltypes, from = summary_filtered_df$pair_cell.types, to = as.vector(summary_filtered_df$Cell_type.source))
celltype_target_vec <- mapvalues(x = interaction_celltypes, from = summary_filtered_df$pair_cell.types, to = as.vector(summary_filtered_df$Cell_type.target))
# celltype_text_vec <- paste0(celltype_source_vec, "->", celltype_target_vec)
rowanno_obj1 <- rowAnnotation(Celltype.Ligand = anno_simple(x = vector(mode = "numeric", length = length(celltype_source_vec)), 
                                                            gp = gpar(fill = NA, col = NA), width = unit(10, "mm"), 
                                                            pch = 21, 
                                                            pt_size = unit(7, "mm"), 
                                                            pt_gp = gpar(fill = colors_celltype[celltype_source_vec], col = NA)), 
                              Gene_source = anno_text(x = gene_source_vec,
                                                      gp = gpar(fontface = "italic", fontsize = 20, border = NA),
                                                      location = 0.5,
                                                      just = "center"),
                              Text_arrow = anno_text(x = rep(x = "->", length(celltype_source_vec)),
                                                     gp = gpar(fontsize = 17, border = NA),
                                                     location = 0.5,
                                                     just = "center"),
                              Celltype.Receptor = anno_simple(x = vector(mode = "numeric", length = length(celltype_target_vec)),
                                                              gp = gpar(fill = NA, col = NA), width = unit(10, "mm"), 
                                                              pch = 21, 
                                                              pt_size = unit(7, "mm"), 
                                                              pt_gp = gpar(fill = colors_celltype[celltype_target_vec], col = NA)),
                              Gene_target = anno_text(x = gene_target_vec,
                                                      gp = gpar(fontface = "italic", fontsize = 20, border = NA),
                                                      location = 0.5,
                                                      just = "center"),
                              annotation_name_side = "bottom")

# make column split -------------------------------------------------------
sampletypes_vec <- mapvalues(x = ids_aliquot, from = specimen_clinical_df$Aliquot.snRNA.WU, to = as.vector(specimen_clinical_df$Sample_Type))
sampletypes_factor <- factor(x = sampletypes_vec, levels = c("Tumor", "Normal"))

# make column annotation --------------------------------------------------
colanno_obj1 <- HeatmapAnnotation(Sample_Type = anno_simple(x = sampletypes_vec, col = colors_sampletype[sampletypes_vec]),
                                  annotation_name_side = "left")

# plot hetamp body --------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat,
                             col = colors_heatmapbody,
                             na_col = color_na, 
                             border = "grey70",
                             ## row
                             row_names_side = "left",
                             row_names_max_width = unit(12, "cm"),
                             left_annotation = rowanno_obj1,
                             show_row_names = F,
                             show_row_dend = F, cluster_rows = T, 
                             row_title_side = "right",
                             ## column
                             show_column_dend = F, cluster_columns = T, show_column_names = F,
                             column_split = sampletypes_factor, column_title = NULL, 
                             top_annotation = colanno_obj1,
                             ## others
                             show_heatmap_legend = F)
p
# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_heatmapbody_legand, 
         title = "Interaction strength", 
         title_gp = gpar(fontsize = 20),
         legend_width = unit(4, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"),
  Legend(labels = names(colors_sampletype),
         title = "Sample type",
         legend_gp = gpar(fill = colors_sampletype), 
         labels_gp = gpar(fontsize=15), title_gp = gpar(fontsize = 20),
         grid_height = unit(0.75, "cm"),
         grid_width = unit(0.75, "cm"), nrow = 1))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "druggable_interactions", ".png")
png(file2write, width = 1500, height = 700, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

file2write <- paste0(dir_out, "druggable_interactions", ".pdf")
pdf(file2write, width = 10, height = 4.5, useDingbats = F)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

