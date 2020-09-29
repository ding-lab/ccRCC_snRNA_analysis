# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 4
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_byindividualcluster_bycellgroup7_byaliquot_on_katmai/20200917.v1/avgexp.SCT.bycellgroup.byaliquot.bycluster.31_aliquot_integration.20200917.v1.tsv", data.table = F)
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20200917.v2/31Aliquot.Barcode2CellType.20200917.v2.tsv", data.table = F)
## barcode 2 individual cluster id
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample/20200717.v1/Barcode2MetaData.20200717.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# specify genes to filter -------------------------------------------------
## input kidney-specific EMT genes
# emt_genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/combine_pt_with_emt_markers/20200911.v1/Kidney_Specific_EMT_Genes.20200911.v1.tsv")
emt_genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/combine_pt_with_emt_markers_all/20200920.v1/Kidney_Specific_EMT_Genes.20200920.v1.tsv")
## add name for the marker groups
emt_genes_df <- emt_genes_df %>%
  mutate(Text_Gene_Group = ifelse(Gene_Group2 == "Tumor cells", 
                                  "Tumor-cell\nmarkers", paste0(Gene_Group2, "\nmarkers")))

genes2filter <- emt_genes_df$Gene
# genes2filter <- emt_genes_df$Gene[emt_genes_df$Gene_Group2 %in% "Mesenchymal"]
# genes2filter <- emt_genes_df$Gene[emt_genes_df$Gene_Group2 %in% c("Mesenchymal", "Epithelial")]
# genes2filter <- emt_genes_df$Gene[emt_genes_df$Gene_Group2 %in% c("Mesenchymal", "Epithelial") & !(emt_genes_df$Gene %in% c("MMP2", "MMP9", "MMP15"))]
# genes2filter <- emt_genes_df$Gene[!(emt_genes_df$Gene %in% c("MMP2", "MMP9", "MMP15"))]
# genes2filter <- emt_genes_df$Gene[emt_genes_df$Gene_Group2 %in% c("Mesenchymal", "Epithelial") & !(emt_genes_df$Gene %in% c("MMP2", "MMP9", "MMP15"))]
# genes2filter <- emt_genes_df$Gene[!(emt_genes_df$Gene %in% c("MMP2", "MMP9", "MMP15"))]

# count cell number and filter clusters -----------------------------------
barcode2celltype_df <- merge(barcode2celltype_df, barcode2cluster_df, by.x = c("orig.ident", "individual_barcode"), by.y = c("aliquot", "individual_barcode"), all.x = T)
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_bycluster_bycellgroup_byaliquot = paste0(orig.ident, "_", seurat_cluster_id, "_",Cell_group7))
cellcount_bycluster_df <- barcode2celltype_df %>%
  select(id_bycluster_bycellgroup_byaliquot) %>%
  table() %>%
  as.data.frame() %>%
  rename(id_bycluster_bycellgroup_byaliquot_original = ".") %>%
  mutate(id_bycluster_bycellgroup_byaliquot = gsub(x = id_bycluster_bycellgroup_byaliquot_original, pattern = " |\\-", replacement = "."))

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycluster_bycellgroup_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  mutate(aliquot = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,1]) %>%
  mutate(id_cluster = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,2]) %>%
  mutate(cellgroup = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,3])
plot_data_long_df$Cell_count <- mapvalues(x = plot_data_long_df$id_bycluster_bycellgroup_byaliquot, from = cellcount_bycluster_df$id_bycluster_bycellgroup_byaliquot, to = as.vector(cellcount_bycluster_df$Freq))
plot_data_long_df$Cell_count <- as.numeric(as.vector(plot_data_long_df$Cell_count))
plot_data_long_df <- plot_data_long_df %>%
  dplyr::filter(Cell_count >= 30) %>%
  dplyr::filter(cellgroup %in% c("Tumor.cells", "Transitional.cells", "Tumor.like.cells"))
plot_data_long_df$id_aliquot_wu <- mapvalues(x = plot_data_long_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
plot_data_long_df <- plot_data_long_df %>%
  dplyr::mutate(id_bycluster_bycellgroup_byaliquot_new = paste0(id_aliquot_wu, "_", id_cluster, "_", cellgroup))
## make matrix
plot_data_wide_df <- dcast(data = plot_data_long_df, formula = V1 ~ id_bycluster_bycellgroup_byaliquot_new, value.var = "value")
plot_data_mat <- as.matrix(plot_data_wide_df[,-1])
rownames(plot_data_mat) <- plot_data_wide_df$V1
plot_data_mat[1:5, 1:5]

# filter genes based on variation -----------------------------------------
sd_bygene_df <- data.frame(SD = apply(plot_data_mat,1, sd, na.rm = TRUE), gene = rownames(plot_data_mat))
sd_bygene_df$Cell_Group2 <- mapvalues(x = sd_bygene_df$gene, from = emt_genes_df$Gene, to = as.vector(emt_genes_df$Gene_Group2))
sd_bygene_df <- sd_bygene_df %>%
  arrange(desc(SD))
genes_plot_mesenchymal <- as.vector(sd_bygene_df$gene[sd_bygene_df$Cell_Group2 == "Mesenchymal"])
genes_plot_mesenchymal
# sd_bygene_other_df <- sd_bygene_df %>%
#   filter(Cell_Group2 != "Mesenchymal") %>%
#   group_by(Cell_Group2) %>%
#   top_n(n = 2, wt = SD)
# genes_plot_other <- as.vector(sd_bygene_other_df$gene)
genes_plot_tumormarkers <- c("PAX8", "PAX2", "CA9")
sd_bygene_epithelial_df <- sd_bygene_df %>%
  filter(Cell_Group2 %in% c("Epithelial", "Proximal tubule")) %>%
  arrange(desc(SD))
genes_plot_epithelial <- head(x = as.vector(sd_bygene_epithelial_df$gene), n = 10)
genes_plot_other <- c(genes_plot_tumormarkers, genes_plot_epithelial)
genes_plot <- c(genes_plot_mesenchymal, genes_plot_other)
plot_data_mat <- plot_data_mat[genes_plot,]
## make mesenchymal score
genes_mesenchymal_score <- c("FN1", "CDH2", "VIM", "FOXC2", "SNAI2")
scores_mesenchymal <- colMeans(plot_data_mat[genes_mesenchymal_score,])
## make mesenchymal score
# genes_epithelial_score <- c("KRT19", "CLDN10", "GPX3", "SLC5A12")
genes_epithelial_score <- head(x = genes_plot_epithelial, n = 5)
scores_epithelial <- colMeans(plot_data_mat[genes_epithelial_score,])

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
# summary(as.vector(plot_data_mat))
colors_heatmapbody = colorRamp2(c(-1.5, 
                                  0, 
                                  1.5), 
                                c(color_blue, "white", color_red))
## make colors for mesenchymal score
summary(scores_mesenchymal)
color_orange <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[5]
color_purple <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[4]
colors_scores_mesenchymal = colorRamp2(c(-0.5, 0, 0.5), 
                                c(color_purple, "white", color_orange))
## make colors for epithelial score
summary(scores_epithelial)
color_yellow <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")[1]
color_blue2 <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")[9]
colors_scores_epithelial = colorRamp2(c(-1, 0, 1), 
                                       c(color_yellow, "white", color_blue2))

# get column ids ----------------------------------------------------------
columnnames_plot <- colnames(plot_data_mat)
ids_aliquot_wu <- str_split_fixed(string = columnnames_plot, pattern = "_", n = 3)[,1]; ids_aliquot_wu
ids_cluster <- str_split_fixed(string = columnnames_plot, pattern = "_", n = 3)[,2]; ids_cluster

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = genes_plot, from = emt_genes_df$Gene, to = as.vector(emt_genes_df$Text_Gene_Group))
row_split_vec
row_split_factor <- factor(x = row_split_vec, levels = c("Mesenchymal\nmarkers", "Epithelial\nmarkers", "Proximal tubule\nmarkers", "Tumor-cell\nmarkers"))

# make column annotation --------------------------------------------------
## make highlighted samples
index_highlight <- which(columnnames_plot %in% c("C3N-01200-T2_5_Transitional.cells", "C3N-01200-T1_5_Tumor.cells", "C3N-01200-T3_2_Transitional.cells", "C3N-00495-T1_10_Transitional.cells",
                                                 "C3L-00079-T1_7_Transitional.cells", "C3L-00790-T1_7_Transitional.cells", "C3N-00495-T1_9_Transitional.cells"))
texts_aliquot_cluster <- paste0(ids_aliquot_wu, "_C", ids_cluster)
texts_highlight <- texts_aliquot_cluster[index_highlight]; texts_highlight
## make column annotation object
colanno_obj = HeatmapAnnotation(MesenchymalScore = anno_simple(x = scores_mesenchymal, col = colors_scores_mesenchymal, height = unit(1, "cm")),
                                EpithelialScore = anno_simple(x = scores_epithelial, col = colors_scores_epithelial, height = unit(1, "cm")),
                                link = anno_mark(at = index_highlight, labels = texts_highlight, labels_gp = gpar(fontsize = 15), side = "bottom"),
                                annotation_name_gp = gpar(fontsize = 20, fontface = "italic"), annotation_name_side = "left")
# make column order --------------------------------------------------
column_order_vec <- order(scores_mesenchymal, decreasing = T)

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na, border = "black",
                             show_row_names = T, row_names_gp = gpar(fontsize = 15, fontface = "bold"),
                             row_split = row_split_factor,
                             row_title_rot = 0, row_title_gp = gpar(fontsize = 20, fontface = "bold"),
                             # row_labels = factor_cellgroup,
                             cluster_row_slices = F, show_row_dend = F, 
                             # column_km = 8, column_km_repeats = 150, 
                             show_column_dend = F, cluster_columns = F, 
                             column_order = column_order_vec,
                             bottom = colanno_obj, show_column_names = F, column_title = NA,
                             show_heatmap_legend = F)
p


# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_heatmapbody, 
         title = "Scaled snRNA expression", 
         title_gp = gpar(fontsize = 15, fontface = "bold"),
         legend_width = unit(6, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_scores_mesenchymal, 
         title = "Mesenchymal score", 
         title_gp = gpar(fontsize = 15, fontface = "bold"),
         legend_width = unit(6, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_scores_epithelial, 
         title = "Epithelial score", 
         title_gp = gpar(fontsize = 15, fontface = "bold"),
         legend_width = unit(6, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "EMT_Genes_by_tumorcluster", ".png")
png(file2write, width = 3000, height = 1200, res = 150)
draw(object = p, 
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()

