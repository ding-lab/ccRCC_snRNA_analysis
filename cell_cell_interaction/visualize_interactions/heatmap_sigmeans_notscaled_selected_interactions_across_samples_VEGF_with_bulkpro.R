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
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/rank_across_celltypes_by_pair_by_avgsigmean/20200925.v3/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
## input summary for cell-cell interactin
# summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/filter_interactions/filter_druggable_related_pair_celltypes/20200924.v1/filtered_druggale_related_pair_celltypes.tsv")
## input cell-cell interaction values by case
# cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/other/scale_pair_cell_types_vs_sample_sig_mean/20200924.v1/cellphonedb.pair_cell.types_vs_sample.sig_mean.scaled_by_sample_bypair..tsv")
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/other/scale_pair_cell_types_vs_sample_sig_mean/20200924.v1/cellphonedb.pair_cell.types_vs_sample.sig_mean.tsv")
## input sample meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input clinical info
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/extract_specimen_clinical_data/20200717.v1/snRNA_ccRCC_Specimen_Clinicl_Data.20200717.v1.tsv")
## input protein data
protein_df <- fread("./Resources/Bulk_Processed_Data/Protein/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)

# specify pairs to filter -------------------------------------------------
## get genes to filter
summary_filtered_df <- summary_df %>%
  filter(pair_cell.types %in% c("VEGFA_FLT1.Tumor cells|Endothelial cells", "VEGFA_KDR.Tumor cells|Endothelial cells", "FLT1_VEGFB.Endothelial cells|Tumor cells")) %>%
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
# plot_data_mat <- plot_data_mat[,!(colnames(plot_data_mat) %in% c("C3L-00359-T1"))]
plot_data_mat[1:3, 1:10]
summary(as.vector(plot_data_mat))
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
## cell type colors
colors_celltype <- RColorBrewer::brewer.pal(n = 4, name = "Set1")
names(colors_celltype) <- c("Tumor cells", "Endothelial cells", "Myofibroblasts", "Fibroblasts")
## specify color for NA values
color_na <- "grey70"
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
## make color for the protein
colors_protein = colorRamp2(c(-1, 
                              0, 
                              1), 
                            c(color_blue, "white", color_red))


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
celltype_text_vec <- paste0(celltype_source_vec, "->", celltype_target_vec)
rowanno_obj1 <- rowAnnotation(Genes = anno_text(x = gene_text_vec, gp = gpar(fontface = "italic", fontsize = 15)),
                              CellTypes = anno_text(x = celltype_text_vec, gp = gpar(fontface = "italic", fontsize = 15)),
                              annotation_name_side = "bottom")

# make column annotation --------------------------------------------------
vegfa_protein_vec <- as.vector(protein_mat2["VEGFA",]); names(vegfa_protein_vec) <- ids_aliquot[!is.na(ids_aliquot_bulk)]; vegfa_protein_vec <- vegfa_protein_vec[ids_aliquot]
vegfb_protein_vec <- as.vector(protein_mat2["VEGFB",]); names(vegfb_protein_vec) <- ids_aliquot[!is.na(ids_aliquot_bulk)]; vegfb_protein_vec <- vegfb_protein_vec[ids_aliquot]
flt1_protein_vec <- as.vector(protein_mat2["FLT1",]); names(flt1_protein_vec) <- ids_aliquot[!is.na(ids_aliquot_bulk)]; flt1_protein_vec <- flt1_protein_vec[ids_aliquot]
kdr_protein_vec <- as.vector(protein_mat2["KDR",]); names(kdr_protein_vec) <- ids_aliquot[!is.na(ids_aliquot_bulk)]; kdr_protein_vec <- kdr_protein_vec[ids_aliquot]

colanno_obj1 <- HeatmapAnnotation(VEGFA.bulk.protein = anno_simple(x = vegfa_protein_vec, col = colors_protein),
                                  FLT1.bulk.protein = anno_simple(x = flt1_protein_vec, col = colors_protein),
                                  VEGFB.bulk.protein = anno_simple(x = vegfb_protein_vec, col = colors_protein),
                                  KDR.bulk.protein = anno_simple(x = kdr_protein_vec, col = colors_protein),
                                  annotation_name_side = "left")

# make column split -------------------------------------------------------
sampletypes_vec <- mapvalues(x = ids_aliquot, from = specimen_clinical_df$Aliquot.snRNA.WU, to = as.vector(specimen_clinical_df$Histologic_Type))

# plot hetamp body --------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat,
                             col = colors_heatmapbody,
                             na_col = color_na, border = "black", 
                             ## row
                             row_names_side = "left",
                             row_names_max_width = unit(12, "cm"),
                             left_annotation = rowanno_obj1,
                             show_row_names = F,
                             show_row_dend = F, cluster_rows = T, 
                             row_title_side = "right",
                             ## column
                             show_column_dend = F, cluster_columns = T, 
                             column_split = sampletypes_vec, column_title = NULL, 
                             bottom_annotation = colanno_obj1,
                             show_heatmap_legend = F)
p
# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_heatmapbody_legand, 
         title = "Scaled interaction strength", 
         title_gp = gpar(fontsize = 10, fontface = "bold"),
         legend_width = unit(4, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_protein, 
         title = "Protein abundance", 
         title_gp = gpar(fontsize = 10, fontface = "bold"),
         legend_width = unit(4, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "druggable_interactions", ".png")
png(file2write, width = 1500, height = 700, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

file2write <- paste0(dir_out, "druggable_interactions", ".pdf")
pdf(file2write, width = 10, height = 3.75, useDingbats = F)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

