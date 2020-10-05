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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/rank_across_celltypes_by_pair_by_avgsigmean/20200925.v3/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
## input summary for cell-cell interactin
# summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/filter_interactions/filter_druggable_related_pair_celltypes/20200924.v1/filtered_druggale_related_pair_celltypes.tsv")
## input cell-cell interaction values by case
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/other/scale_pair_cell_types_vs_sample_sig_mean/20200924.v1/cellphonedb.pair_cell.types_vs_sample.sig_mean.scaled_by_sample_bypair..tsv")
## input sample meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input clinical info
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/extract_specimen_clinical_data/20200717.v1/snRNA_ccRCC_Specimen_Clinicl_Data.20200717.v1.tsv")

# specify pairs to filter -------------------------------------------------
## get genes to filter
paired_cellgroups.general_process <- "Tumor&Stroma"
summary_filtered_df1 <- summary_df %>%
  filter(paired_cellgroups.general == paired_cellgroups.general_process) %>%
  filter(rank_byinteraction_acrosscelltypes == 1) %>%
  filter(rank_genesource_acrosscelltypes == 1 & rank_genetarget_acrosscelltypes == 1) %>%
  filter(!(interacting_pair %in% c("TNC_aVb6 complex", "PVR_NECTIN3", "FGF7_FGFR4"))) %>%
  mutate(paired_celltypes_group = paste0("Tumor&",ifelse(Cell_type.source == "Tumor cells", Cell_type.target, Cell_type.source))) %>%
  arrange(desc(avg_sig_mean), number_sig_cases)
paired_cellgroups.general_process <- "Stroma&Stroma"
summary_filtered_df2 <- summary_df %>%
  filter(paired_cellgroups.general == paired_cellgroups.general_process) %>%
  filter(rank_byinteraction_acrosscelltypes == 1) %>%
  filter(rank_genesource_acrosscelltypes == 1 & rank_genetarget_acrosscelltypes == 1) %>%
  filter(!is.na(gene.source.druggable) | !is.na(gene.target.druggable)) %>%
  filter(!is_integrin) %>%
  filter(!(interacting_pair %in% c("PlexinA2_complex1_SEMA3A", "PlexinA1_complex1_SEMA6D", "ACVR_1B2A receptor_INHBA", "EPOR_KITLG", "FGFR1_FGF7", "FGF1_FGFR1"))) %>%
  mutate(paired_celltypes_group = ifelse(Cell_type.target == "Endothelial cells", paste0("Endothelial cells", "&", Cell_type.source), paste0(Cell_type.source, "&", Cell_type.target))) %>%
  arrange(desc(avg_sig_mean), number_sig_cases)
summary_filtered_df <- rbind(summary_filtered_df1, summary_filtered_df2)
summary_filtered_df <- summary_filtered_df %>%
  filter(interacting_pair %in% c("VEGFA_FLT1", "VEGFA_KDR", "FLT4_VEGFC", "ANGPT2_TEK", "TEK_ANGPT1", "PDGFRB_PDGFD")) %>%
  arrange(desc(avg_sig_mean))
interacting_pairs_process <- summary_filtered_df$interacting_pair
interacting_pairs_process
length(interacting_pairs_process)
paired_celltypes_process <- summary_filtered_df$paired_celltypes

# specify the aliquot ids (ordered) ---------------------------------------

# make data matrix for the heatmap body -------------------------------------------
plot_data_df <- cellphone_df %>%
  filter(V1 %in% summary_filtered_df$pair_cell.types)
plot_data_mat <- as.matrix(plot_data_df[, -1])
rownames(plot_data_mat) <- plot_data_df$V1
## remove the sample without endothelial cells
plot_data_mat <- plot_data_mat[,!(colnames(plot_data_mat) %in% c("C3L-00359-T1"))]
plot_data_mat[1:5, 1:10]
summary(as.vector(plot_data_mat))
##
plot_data_mat[is.na(plot_data_mat)] <- -3
ids_aliquot <- colnames(plot_data_mat)
interaction_celltypes <- rownames(plot_data_mat)

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
colors_heatmapbody = colorRamp2(c(-3,
                                  -2, 
                                  0, 
                                  2), 
                                c(color_na, color_blue, "white", color_red))
colors_heatmapbody_legand <- colorRamp2(c(-2, 
                                          0, 
                                          2), 
                                        c(color_blue, "white", color_red))

## make colors for average sig means
colors_sig_means <- colorRamp2(c(0, 2, 4, 6, 8, 10), 
                               brewer.pal(n = 6, name = "YlGn"))

# make row annotation -----------------------------------------------------
avg_sig_mean_vec <- mapvalues(x = interaction_celltypes, from = summary_df$pair_cell.types, to = as.vector(summary_df$avg_sig_mean))
avg_sig_mean_vec <- as.numeric(avg_sig_mean_vec)
avg_sig_mean_vec
summary(avg_sig_mean_vec)
## make text
gene_source_vec <- mapvalues(x = interaction_celltypes, from = summary_filtered_df$pair_cell.types, to = as.vector(summary_filtered_df$gene.source))
gene_target_vec <- mapvalues(x = interaction_celltypes, from = summary_filtered_df$pair_cell.types, to = as.vector(summary_filtered_df$gene.target))
## make cell types
celltype_source_vec <- mapvalues(x = interaction_celltypes, from = summary_filtered_df$pair_cell.types, to = as.vector(summary_filtered_df$Cell_type.source))
celltype_target_vec <- mapvalues(x = interaction_celltypes, from = summary_filtered_df$pair_cell.types, to = as.vector(summary_filtered_df$Cell_type.target))
rowanno_obj <- rowAnnotation(Gene_source = anno_text(x = gene_source_vec, gp = gpar(fill = colors_celltype[celltype_source_vec])),
                             Gene_target = anno_text(x = gene_target_vec, gp = gpar(fill = colors_celltype[celltype_target_vec])))

rowanno_obj1 <- rowAnnotation(Celltype_ligand = anno_simple(x = vector(mode = "numeric", length = length(celltype_source_vec)), 
                                                            gp = gpar(fill = "white"),
                                                            pch = 21, pt_gp = gpar(fill = colors_celltype[celltype_source_vec])),
                              Gene_source = anno_text(x = gene_source_vec, gp = gpar(fontface = "italic", fontsize = 15)),
                              Celltype_receptor = anno_simple(x = vector(mode = "numeric", length = length(celltype_target_vec)), 
                                                              gp = gpar(fill = "white"),
                                                              pch = 21, pt_gp = gpar(fill = colors_celltype[celltype_target_vec])),
                              Gene_target = anno_text(x = gene_target_vec, gp = gpar(fontface = "italic", fontsize = 15)),
                              # Pathway = anno_text(x = pathway_vec, gp = gpar(fill = colors_pathway[pathway_vec])),
                              annotation_name_side = "top")

# make row split ----------------------------------------------------------
therapy_category_vec <- mapvalues(x = interaction_celltypes, from = summary_df$pair_cell.types, to = as.vector(summary_df$therapy_category))
therapy_category_vec

# make column annotation --------------------------------------------------


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
  # Legend(col_fun = colors_sig_means, 
  #        title = "Interaction strength score", 
  #        title_gp = gpar(fontsize = 10, fontface = "bold"),
  #        legend_width = unit(4, "cm"),
  #        legend_height = unit(3, "cm"),
  #        direction = "horizontal"),
  Legend(labels = names(colors_celltype),
         title = "Cell_type",
         legend_gp = gpar(fill = colors_celltype)))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "druggable_interactions", ".png")
png(file2write, width = 1000, height = 700, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

file2write <- paste0(dir_out, "druggable_interactions", ".pdf")
pdf(file2write, width = 6.2, height = 5)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

