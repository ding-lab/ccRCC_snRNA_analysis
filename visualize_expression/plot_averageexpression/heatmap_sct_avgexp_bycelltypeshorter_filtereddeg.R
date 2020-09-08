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
genes_df <- fread(input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/filter_markers/filter_markers_wilcox_bygroup/20200904.v1/findallmarkers_wilcox_bycellgroup.pos.logfcthreshold0.1.minpct0.1.mindiffpct0.1.top50avg_logFC.tsv", data.table = F)
## input the average expression calculated (RNA)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_bycelltypeshorter_on_katmai/20200904.v1/averageexpression_SCT_bycelltype.shorter.31_aliquot_integration.20200904.v1.tsv", data.table = F)
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200904.v1/31Aliquot.Barcode2CellType.20200904.v1.tsv", data.table = F)
## input gene annotation
genes_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_ccrcc_pathogenic_pathway_genes/20200907.v1/ccRCC_Pathogenic_Pathways_Genes.20200907.v1.tsv")
## input not scaled average expression
avgexp_baseline_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/format_avgexp_sct_usedata_bycelltypeshorter/20200907.v2/formated.averageexpression.SCT.slotdata.bycelltype.shorter.31_aliquot_integration.20200907.v2.tsv")
## input gene to cell type enriched TFs
gene2motifenrichedtf_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/annotate_markers/annotate_degs_to_celltypespecific_tfs/20200908.v1/CellTypeDEG2CellTypeMotifEnrichedTFs.Annotation.20200908.v1.tsv")

# specify the genes to show -----------------------------------------------
genes2filter <- genes_df$row_name
# genes2filter <- c(genes2filter, ccRCC_drivers)

# calculate the cell counts and specify the cell types to show ------------
cellcount_df <- barcode2celltype_df %>%
  select(Cell_type.shorter) %>%
  table() %>%
  as.data.frame() %>%
  rename(Cell_type.shorter = '.') %>%
  mutate(Keep = (Freq >= 50 & !(Cell_type.shorter %in% c("Unknown", "Tumor-like cells")))) %>%
  # mutate(Keep = (Freq >= 50 & Cell_type.shorter != "Unknown")) %>%
  mutate(colname_celltype = gsub(x = Cell_type.shorter, pattern = " |\\+|\\/|\\-", replacement = "."))
data_col_names.keep <- cellcount_df$colname_celltype[cellcount_df$Keep]
data_col_names.keep
celltype2cellgroup_df <- barcode2celltype_df %>%
  select(Cell_type.shorter, Cell_group.shorter, Cell_group.detailed) %>%
  unique() %>%
  mutate(colname_celltype = gsub(x = Cell_type.shorter, pattern = " |\\+|\\/|\\-", replacement = "."))

# format the column names to only aliquot id ------------------------------
## filtr the rows
plot_data_df <- avgexp_df %>%
  rename(gene = V1) %>%
  filter(gene %in% genes2filter)
## remove teh prefix from the column names
data_col_names <- colnames(plot_data_df)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
data_col_names.changed
## rename the data frame
colnames(plot_data_df) <- c("gene", data_col_names.changed)
## filter the columns and make data matrix
plot_data_mat <- t(as.matrix(plot_data_df[, data_col_names.keep]))
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
factor_cellgroup <- factor(x = vec_cellgroup, levels = c("Nephron_Epithelium", "Stroma", "Immune"))

# make annotation for the genes -------------------------------------------
## annotate the pathways
vec_is_vhlhif <- as.character(genes_plot %in% genes_anno_df$target_genesymbol[genes_anno_df$pathway_name == "VHL-HIF"])
vec_is_epigentic_smg_related <- as.character(genes_plot %in% genes_anno_df$target_genesymbol[genes_anno_df$pathway_name == "Epigenetic machinary"])
vec_is_pi3kmtor <- as.character(genes_plot %in% genes_anno_df$target_genesymbol[genes_anno_df$pathway_name == "PI3K-AKT-mTOR Pathway"])
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
idx_is_highlighted <- which((vec_is_vhlhif == "TRUE" | vec_is_epigentic_smg_related == "TRUE" | vec_is_pi3kmtor == "TRUE" | vec_targetof_motifenrichedtfs == "Stimulated Target of Cell-Type-Motif-Enriched TFs"))
genes_highlighted <- genes_plot[idx_is_highlighted]
## make annotation object
colanno_obj = HeatmapAnnotation(Baseline_Expression = anno_simple(x = vec_avgexp_baseline, col = col_baselineexp), 
                                VHL_HIF_Pathway = anno_simple(x = vec_is_vhlhif, col = colors_vhlhif),
                                Epigenetic_machinary_related = anno_simple(x = vec_is_epigentic_smg_related, col = colors_vhlhif),
                                PI3K_mTOR_Pathway = anno_simple(x = vec_is_pi3kmtor, col = colors_vhlhif),
                                Is_Targetof_MotifEnrichedTFs_byCellGroup = anno_simple(x = vec_targetof_motifenrichedtfs, col = colors_targetof_motifenrichedtfs),
                                link = anno_mark(at = idx_is_highlighted, labels = genes_highlighted, labels_gp = gpar(fontsize = 10), side = "bottom"),
                                annotation_name_gp = gpar(fontsize = 10))

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
             row_split = factor_cellgroup, cluster_row_slices = F, show_row_dend = F, row_title_rot = 0, row_title_gp = gpar(fontsize = 10),
             show_column_names = T, column_names_gp = gpar(fontsize = 5), column_km = 6, column_km_repeats = 100, show_column_dend = F, 
             bottom_annotation = colanno_obj, 
             show_heatmap_legend = F)
p
## save heatmap as png
png(filename = paste0(dir_out, "heatmap", ".png"), 
    width = 2000, height = 1000, res = 150)
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




