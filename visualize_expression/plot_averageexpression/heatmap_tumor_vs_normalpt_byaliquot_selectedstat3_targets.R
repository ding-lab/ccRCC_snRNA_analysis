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
avgexp_bycelltype_byaliquot_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_bycelltypedetailed_byaliquot_on_katmai/20200413.v1/averageexpression_bycelltypedetailed.30_aliquot_integration.20200413.v1.tsv", data.table = F)
## input TF integrations
tf_tab <- fread(input = "./Resources/PPI/TF_interactions.txt", data.table = F)
## input barcode 2 cell type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200410.v1/30_aliquot_integration.barcode2celltype.20200410.v1.tsv", data.table = F)
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200413.v1/meta_data.20200413.v1.tsv", data.table = F)
## load bulk protein data
protein_df <- fread("./Resources/Bulk_Processed_Data/Protein/CCRCC_PRO_tumor_PGDAC_MD_MAD_partID.txt", data.table = F)
## input deg summary file
deg_summary_df <- fread("./Resources/Analysis_Results/findmarkers/tumor_vs_normal/examine_degs_roc_for_tumor_vs_normalpt/20200421.v1/summary_findmarkers_roc.tumor_vs_normalpt.20200421.v1.tsv", data.table = F)

# set plot gene list ------------------------------------------------------
genes_plot <- tf_tab$target_genesymbol[tf_tab$source_genesymbol == "STAT3" & tf_tab$is_stimulation == 1]
genes_plot
# genes_plot <- c(genes_plot, "RORA", "VEGFA", "HIF1A", "NR3C1", "STAT3", "PFKFB3", "TFRC", "EGLN3", "POU5F1")
genes_plot <- c("HIF1A", "STAT3", genes_plot)
genes_plot <- intersect(genes_plot, deg_summary_df$gene_symbol)
genes_plot

# count number of cells per celltype per aliquot --------------------------
count_bycelltype_byaliquot_df <- data.frame(table(barcode2celltype_df[, c("Cell_type.detailed", "orig.ident")]))
unique(count_bycelltype_byaliquot_df$Cell_type.detailed)
count_bycelltype_byaliquot_df$celltype <- gsub(x = count_bycelltype_byaliquot_df$Cell_type.detailed, pattern = '[/ +]', replacement = ".")
count_bycelltype_byaliquot_df$celltype <- gsub(x = count_bycelltype_byaliquot_df$celltype, pattern = '\\-', replacement = ".")
count_bycelltype_byaliquot_df <- count_bycelltype_byaliquot_df %>%
  mutate(id_aliquot_celltype = paste0(orig.ident, "_", celltype))
count_bycelltype_byaliquot_filtered <- count_bycelltype_byaliquot_df %>%
  filter(celltype %in% c("Tumor.cells", "Proximal.tubule")) %>%
  filter(Freq >= 10)

# make matrix for heatmap body --------------------------------------------
plot_data_df <- avgexp_bycelltype_byaliquot_df
## format the column names to only aliquot id + cell type
data_col_names <- colnames(plot_data_df)[-1]
data_col_names
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
data_col_names.changed
## rename the data frame
colnames(plot_data_df) <- c("gene", data_col_names.changed)
## get the column names for tumor cells and normal epithelial cells with enough cells
data_col_names.keep <- count_bycelltype_byaliquot_filtered$id_aliquot_celltype
# data_col_names.keep <- data_col_names[grepl(pattern = "Normal.epithelial.cells", x = data_col_names) | grepl(pattern = "Tumor.cells", x = data_col_names)]
data_col_names.keep
## get the gene names for HIF targets
data_row_names <- plot_data_df$gene
data_row_names.keep <- intersect(genes_plot, data_row_names)
data_row_names.keep
## filter cell group down to epithelial cells
## reformat data frame to matrix
plot_data_mat <- as.matrix(plot_data_df[,data_col_names.keep])
plot_data_mat %>% head()
## add row names
rownames(plot_data_mat) <- data_row_names
plot_data_mat %>% head()
## filter down to HIF target genes
plot_data_mat <- plot_data_mat[data_row_names.keep,]
### get aliquot ids and case ids
ids_aliquot_celltype <- data_col_names.keep
ids_aliquot <- str_split_fixed(string = ids_aliquot_celltype, pattern = "_", n = 2)[,1]
ids_aliquot_wu <- mapvalues(x = ids_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
ids_celltype <- str_split_fixed(string = ids_aliquot_celltype, pattern = "_", n = 2)[,2]
ids_case <- mapvalues(x = ids_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
# case_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

# make top annotation -----------------------------------------------------
## make colors for cell types
colors_celltype <- c("green", "yellow")
names(colors_celltype) <- c("Tumor.cells", "Proximal.tubule")
## make protein expression of RORA and HIF1A
protein_selectedgenes_df <- protein_df %>%
  filter(Gene %in% c("RORA", "HIF1A"))
protein_selectedgenes_mat <- protein_selectedgenes_df[,ids_case]
rownames(protein_selectedgenes_mat) <- protein_selectedgenes_df$Gene
colnames(protein_selectedgenes_mat) <- ids_aliquot_wu
protein_selectedgenes_mat[,ids_celltype == "Proximal.tubule"] <- NA
## make top annotation
top_col_anno <- HeatmapAnnotation(Cell.Type = anno_simple(x = ids_celltype, col = colors_celltype))

# make row annotation -----------------------------------------------------
is_hif1a_target <- data_row_names.keep %in% hiftargets_df$target_genesymbol[hiftargets_df$source_genesymbol == "HIF1A"]
is_hif1a_target
colors_is_hif1a_target <- c("black", "white")
names(colors_is_hif1a_target) <- c("TRUE", "FALSE")

is_epas1_target <- data_row_names.keep %in% hiftargets_df$target_genesymbol[hiftargets_df$source_genesymbol == "EPAS1"]
is_epas1_target
colors_is_epas1_target <- c("black", "white")
names(colors_is_epas1_target) <- c("TRUE", "FALSE")
row_anno = rowAnnotation(Is_HIF1A_Target = anno_simple(x = as.character(is_hif1a_target), col = colors_is_hif1a_target),
                         Is_EPAS1_Target = anno_simple(x = as.character(is_epas1_target), col = colors_is_epas1_target))

# make heatmap body ------------------------------------------------------------
heatmapbody_color_fun <- colorRamp2(c(quantile(plot_data_mat, 0.05, na.rm=T), 
                                      quantile(plot_data_mat, 0.5, na.rm=T), 
                                      quantile(plot_data_mat, 0.95, na.rm=T)),
                                    c("blue", "white", "red"))

p <- Heatmap(matrix = plot_data_mat,
             col = heatmapbody_color_fun, 
             right_annotation = row_anno, cluster_rows = T,
             # bottom_annotation = bottom_col_anno,
             top_annotation = top_col_anno,
             column_labels = ids_aliquot_wu,
             show_heatmap_legend = F)
p
# make legend -------------------------------------------------------------
## make legend for top annotation
annotation_lgd = list(
  Legend(col_fun = heatmapbody_color_fun, 
         title = "Average expression value\nby cell type (Normalized)", 
         direction = "vertical"),
  Legend(labels = names(colors_celltype), 
         title = "Cell type", 
         legend_gp = gpar(fill = colors_celltype)))

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "selected_stat3_targets.", "expression.", "tumor_vs_normal_proximal_tubule.", run_id, ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()
