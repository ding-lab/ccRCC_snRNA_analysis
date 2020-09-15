# Yige Wu @WashU Aug 2020

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
## input the DEG  list
# genes_df <- fread(input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/filter_markers/filter_markers_wilcox_bygroup/20200904.v1/findallmarkers_wilcox_bycellgroup.pos.logfcthreshold0.1.minpct0.1.mindiffpct0.1.top50avg_logFC.tsv", data.table = F)
genes_df <- fread(input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/filter_markers/filter_markers_wilcox_bygroup/20200908.v1/findallmarkers_wilcox_bycellgroup.pos.logfcthreshold0.1.minpct0.1.mindiffpct0.1.Top50avg_logFC.tsv", data.table = F)
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/format_expression/combine_avexp_sct_usescale_immunegroups_with_stroma_and_epithelium/20200908.v1/Combined.Average_Expression.SCT.UseScale.20200908.v1.tsv", data.table = F)
## input cell type 2 cell group table
celltype2cellgroup_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/format_expression/combine_avexp_sct_usescale_immunegroups_with_stroma_and_epithelium/20200908.v1/Combined.CellTypes2CellGroup.20200908.v1.tsv")
## input gene annotation
genes2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_ccrcc_pathogenic_pathway_genes/20200908.v1/ccRCC_Pathogenic_Pathways_Genes.20200908.v1.tsv")
## input not scaled average expression
avgexp_baseline_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/format_expression/format_avgexp_sct_usedata_bycelltypeshorter/20200907.v2/formated.averageexpression.SCT.slotdata.bycelltype.shorter.31_aliquot_integration.20200907.v2.tsv")
## input gene to cell type enriched TFs
gene2motifenrichedtf_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/annotate_markers/annotate_degs_to_celltypespecific_tfs/20200908.v2/CellTypeDEGTop502CellTypeMotifEnrichedTFs.Annotation.20200908.v2.tsv")
gene2motifenrichedtf_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/annotate_markers/annotate_degs_to_celltypespecific_tfs/20200908.v2/CellTypeDEGTop502CellTypeMotifEnrichedTFs.20200908.v2.tsv")

# specify the genes to show -----------------------------------------------
genes2filter <- genes_df$row_name
# genes2filter <- c(genes2filter, ccRCC_drivers)

# format the column names to only aliquot id ------------------------------
## filtr the rows
plot_data_df <- avgexp_df %>%
  filter(gene %in% genes2filter)
## filter the columns and make data matrix
plot_data_mat <- t(as.matrix(plot_data_df[,-1]))
## add gene name
colnames(plot_data_mat) <- plot_data_df$gene
plot_data_mat %>% head()

# get dimension names -----------------------------------------------------
genes_plot <- colnames(plot_data_mat)
genes_plot
celltypes_plot <- rownames(plot_data_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color for the pathway annotation
colors_vhlhif <- c("TRUE" = "black", "FALSE" = "white")
## make color function for baseline expression
col_baselineexp <- colorRamp2(c(0, 2), c("white", "orange"))
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
col_fun = colorRamp2(c(-1.5, 
                       0, 
                       1.5), 
                     c(color_blue, "white", color_red))
## make color for the pathway annotation
# RColorBrewer::display.brewer.pal(n = 8, name = "Set1")
colors_targetof_motifenrichedtfs <- c("Stimulated Target of Cell-Type-Motif-Enriched TFs" = RColorBrewer::brewer.pal(n = 9, name = "PuRd")[7], 
                                      "Target of Cell-Type-Motif-Enriched TFs; Direction unknown" = RColorBrewer::brewer.pal(n = 9, name = "PuRd")[3],
                                      "Not a Target" = "white")
colors_targetof_motifenrichedtfs

# make row split for cell types ------------------------------------------
vec_cellgroup <- mapvalues(x = celltypes_plot, from = celltype2cellgroup_df$colname_celltype, to = as.vector(celltype2cellgroup_df$Cell_group.shorter))
vec_cellgroup
vec_cellgroup[vec_cellgroup == "Nephron_Epithelium"] <- "Nephron\nEpithelium"
factor_cellgroup <- factor(x = vec_cellgroup, levels = c("Nephron\nEpithelium", "Stroma", "Immune"))

# make row labels for cell types ------------------------------------------
vec_celltype_labels <- mapvalues(x = celltypes_plot, from = celltype2cellgroup_df$colname_celltype, to = as.vector(celltype2cellgroup_df$Cell_type))
vec_celltype_labels

# make annotation for the genes -------------------------------------------
## annotate the pathways
vec_is_vhlhif <- as.character(genes_plot %in% genes2pathway_df$target_genesymbol[genes2pathway_df$pathway_name == "VHL-HIF"])
vec_is_epigentic_smg_related <- as.character(genes_plot %in% genes2pathway_df$target_genesymbol[genes2pathway_df$pathway_name == "Epigenetic machinary"])
vec_is_pi3kmtor <- as.character(genes_plot %in% genes2pathway_df$target_genesymbol[genes2pathway_df$pathway_name == "PI3K-AKT-mTOR Signaling"])
vec_is_adhesion <- as.character(genes_plot %in% genes2pathway_df$target_genesymbol[genes2pathway_df$pathway_name == "Focal adhesion"])
## annotate the baseline expression
avgexp_baseline_df <- avgexp_baseline_df %>%
  filter(gene %in% genes2filter)
vec_avgexp_baseline <- colMeans(t(as.matrix(avgexp_baseline_df[, celltypes_plot])))
names(vec_avgexp_baseline) <- avgexp_baseline_df$gene
vec_avgexp_baseline <- vec_avgexp_baseline[genes_plot]
summary(vec_avgexp_baseline)
## annotate the genes that are targets of the motif-enriched TFs
vec_targetof_motifenrichedtfs <- gene2motifenrichedtf_anno_df$TF_evidence_level
names(vec_targetof_motifenrichedtfs) <- gene2motifenrichedtf_anno_df$target_genesymbol
vec_targetof_motifenrichedtfs <- vec_targetof_motifenrichedtfs[genes_plot]
vec_targetof_motifenrichedtfs
## annotate genes to highlight
idx_is_highlighted <- which((vec_is_vhlhif == "TRUE" | vec_is_epigentic_smg_related == "TRUE" | vec_is_pi3kmtor == "TRUE" | vec_is_adhesion == "TRUE" | vec_targetof_motifenrichedtfs == "Stimulated Target of Cell-Type-Motif-Enriched TFs"))
genes_highlighted <- genes_plot[idx_is_highlighted]
## make annotation object
colanno_obj = HeatmapAnnotation(#Baseline_Expression = anno_simple(x = vec_avgexp_baseline, col = col_baselineexp), 
                                VHL_HIF_Pathway = anno_simple(x = vec_is_vhlhif, col = colors_vhlhif),
                                Epigenetic_machinary_related = anno_simple(x = vec_is_epigentic_smg_related, col = colors_vhlhif),
                                PI3K_mTOR_Signaling = anno_simple(x = vec_is_pi3kmtor, col = colors_vhlhif),
                                Focal_Adhesion = anno_simple(x = vec_is_adhesion, col = colors_vhlhif),
                                Is_Targetof_MotifEnrichedTFs = anno_simple(x = vec_targetof_motifenrichedtfs, col = colors_targetof_motifenrichedtfs),
                                link = anno_mark(at = idx_is_highlighted, labels = genes_highlighted, labels_gp = gpar(fontsize = 15), side = "bottom"),
                                annotation_name_gp = gpar(fontsize = 15, fontface = "italic"))

# Heatmap -----------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
             # column_labels = ids_aliquot_wu,
             # row_labels = ids_aliquot_wu,
             # right_annotation = row_anno,
             # show_row_names = F, show_column_names = F, 
             # row_split = row_split_factor, cluster_row_slices = F, row_order = row_order_vec,
             # show_row_dend = F, row_title_rot = 0, row_title_side = "right", row_title_gp = gpar(fontsize = 25, fontface = "bold"),
             # row_gap = unit(0, "mm"),
             # column_split = column_split_factor, cluster_column_slices = F, column_order = column_order_vec,
             # show_column_dend = F, column_title_side = "top", column_title_rot = 90, column_title_gp = gpar(fontsize = 25, fontface = "bold"),
             col = col_fun, na_col = color_na,
             show_row_names = T, 
             row_split = factor_cellgroup, 
             row_labels = vec_celltype_labels, row_names_gp = gpar(fontsize = 15, fontface = "bold"),
             cluster_row_slices = F, show_row_dend = F, row_title_rot = 0, row_title_gp = gpar(fontsize = 20),
             show_column_names = F, column_names_gp = gpar(fontsize = 5), column_km = 6, column_km_repeats = 100, show_column_dend = F, 
             bottom_annotation = colanno_obj, 
             show_heatmap_legend = F)
p
## save heatmap as png
png(filename = paste0(dir_out, "heatmap", ".png"), 
    width = 2000, height = 1000, res = 150)
draw(object = p)
dev.off()
## save heatmap as pdf
pdf(paste0(dir_out, "heatmap", ".pdf"), 
    width = 20, height = 6, useDingbats = F)
draw(object = p)
dev.off()

 # legend ------------------------------------------------------------------
list_lgd = list(
  Legend(col_fun = col_fun, 
         title = "Scaled (across all cells) expression values\naveraged within each cell type", 
         legend_width = unit(6, "cm"),
         direction = "horizontal"),
  Legend(col_fun = col_baselineexp, 
         title = "Baseline (not scaled) expression values\naveraged across all cells", 
         legend_width = unit(6, "cm"),
         direction = "horizontal"),
  Legend(title = "Target of cell-type-motif-enriched\ntranscription factor(s)",
         labels = names(colors_targetof_motifenrichedtfs),
         legend_gp = gpar(fill = colors_targetof_motifenrichedtfs)),
  Legend(title = "Member of the ccRCC\nPathogenic Pathways",
         labels = names(colors_vhlhif),
         legend_gp = gpar(fill = colors_vhlhif)))
## create a directory for legends
dir_out_legend <- paste0(dir_out, "Legends/")
dir.create(dir_out_legend)
for (i in 1:length(list_lgd)) {
  file2write <- paste0(dir_out_legend, "Legend", i, ".pdf")
  pdf(file2write,
      width = 4, height = 3)
  draw(list_lgd[[i]])
  dev.off()
  
  file2write <- paste0(dir_out_legend, "Legend", i, ".png")
  png(file2write,
      width = 800, height = 600, res = 150)
  draw(list_lgd[[i]])
  dev.off()
}




