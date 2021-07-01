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
## input emt score pre-calculated
emt_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_emt_scores_tumorcells_per_manual_tumorcluster_scaled/20201201.v1/EMT_scores_by_manual_tumorcluster.20201201.v1.tsv")
## input cell number per cluster
cellnumber_percluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/count_cellnumber_per_manual_cluster_rm_doublet/20210413.v1/CellNumberPerTumorManualCluster.20210413.v1.tsv")
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210624.v1/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

# preprocess-------------------------------------------------
## specify genes to filter 
emt_genes_df <- data.frame(gene = c("VIM", "FN1", "CDH2", "SERPINE1", "TGFBI",
                                    "CUBN", "LRP2", "SLC17A3", "GATM", "GLYAT"),
                           Text_Gene_Group = c(rep("Mesenchymal\nmarkers", 5),
                                               rep("Epithelial/\nproximal-tubule\nmarkers", 5)))
## add name for the marker groups
genes2filter <- emt_genes_df$gene
## identify clusters with sufficient cell number
cluster_pass_df <- cellnumber_percluster_df %>%
  filter(Freq >= 50)%>%
  mutate(colname_exp = gsub(x = id_cluster_uniq,pattern = "\\-", replacement = "."))

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycluster_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  dplyr::filter(id_bycluster_byaliquot %in% cluster_pass_df$colname_exp) %>%
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
colors_heatmapbody = colorRamp2(c(-1, 
                                  0, 
                                  1), 
                                c(color_blue, "white", color_red))
## make colors for mesenchymal score
# summary(scores_mesenchymal)
color_orange <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[5]
color_purple <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[4]
colors_scores_mesenchymal = colorRamp2(c(-0.5, 0, 0.5), 
                                c(color_purple, "white", color_orange))
## make colors for epithelial score
# summary(scores_epithelial)
color_brown <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[7]
color_green <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[3]
colors_scores_epithelial = colorRamp2(c(-1, 0, 1), 
                                       c(color_green, "white", color_brown))
## colors for EMT group
colors_emtgroups <- c("black", "white")
names(colors_emtgroups) <- c("High", "Low")
colors_emtenriched <- c("black", "white smoke")
names(colors_emtenriched) <- c("TRUE", "FALSE")
## make colors for EMT score
colors_emtscores <- colorRamp2(c(-100, 0, 100), 
                               c(color_purple, "white", color_orange))
# get row/column ids ----------------------------------------------------------
columnnames_plot <- colnames(plot_data_mat)
ids_aliquot_wu <- str_split_fixed(string = columnnames_plot, pattern = "_", n = 2)[,1]; ids_aliquot_wu
ids_cluster <- str_split_fixed(string = columnnames_plot, pattern = "_", n = 2)[,2]; ids_cluster
rownames_plot <- rownames(plot_data_mat)

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = rownames_plot, from = emt_genes_df$gene, to = as.vector(emt_genes_df$Text_Gene_Group))
row_split_vec
row_split_factor <- factor(x = row_split_vec, levels = c("Mesenchymal\nmarkers", "Epithelial/\nproximal-tubule\nmarkers", "Tumor-cell\nmarkers"))

# make column annotation --------------------------------------------------
## make cutoff for the EMT-potential-high
emt_scores_df <- emt_scores_df %>%
  mutate(easyid_cluster = paste0(Sample_id, "_", Cluster_id))
emt_group_vec <- mapvalues(x = columnnames_plot, from = emt_scores_df$easyid_cluster, to = as.vector(emt_scores_df$EMT_potential))
scores_mesenchymal <- mapvalues(x = columnnames_plot, from = emt_scores_df$easyid_cluster, to = as.vector(emt_scores_df$Score_mesenchymal)); scores_mesenchymal <- as.numeric(scores_mesenchymal)
scores_epithelial <- mapvalues(x = columnnames_plot, from = emt_scores_df$easyid_cluster, to = as.vector(emt_scores_df$Score_epithelial)); scores_epithelial <- as.numeric(scores_epithelial)
## make 
emt_enriched_vec <- as.character(columnnames_plot %in% enrich_df$cluster_name[enrich_df$EMT])
emt_scores_vec <- mapvalues(x = columnnames_plot, from = enrich_df$cluster_name, to = as.vector(enrich_df$EPITHELIAL_MESENCHYMAL_TRANSITION_Score))
emt_scores_vec <- as.numeric(emt_scores_vec)

## make highlighted samples
index_highlight <- which(columnnames_plot %in% c("C3L.00079.T1_C4", "C3L.00079.T1_C1",  "C3L.00079.T1_C2",  "C3L.00079.T1_C3",  "C3L.00079.T1_C4", 
                                                 "C3N.01200.T2_C1", "C3N.01200.T2_C2",
                                                 "C3N.00242.T1_C1"))
texts_highlight <- columnnames_plot[index_highlight];
## make column annotation object
colanno_obj = HeatmapAnnotation(link = anno_mark(at = index_highlight, labels = texts_highlight, labels_gp = gpar(fontsize = 15), side = "top"),
                                EMTEnriched = anno_simple(x = emt_enriched_vec, col = colors_emtenriched[emt_enriched_vec], height = unit(1, "cm")),
                                EMTScore = anno_simple(x = emt_scores_vec, col = colors_emtscores, height = unit(0.75, "cm")),
                                # EMTPotential = anno_simple(x = emt_group_vec, col = colors_emtgroups[emt_group_vec], height = unit(1, "cm")),
                                # MesenchymalScore = anno_simple(x = scores_mesenchymal, col = colors_scores_mesenchymal, height = unit(0.75, "cm")),
                                # EpithelialScore = anno_simple(x = scores_epithelial, col = colors_scores_epithelial, height = unit(0.75, "cm")),
                                annotation_name_gp = gpar(fontsize = 20, fontface = "bold"), annotation_name_side = "left")

# make column order --------------------------------------------------
emt_group_factor <- factor(x = emt_enriched_vec, levels = c("FALSE", "TRUE"))
column_order_vec <- order(emt_scores_vec, decreasing = T)

# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(title = "EMT gene expression enriched", 
         title_gp = gpar(fontsize = 20),
         labels = names(colors_emtenriched), legend_gp = gpar(fill = colors_emtenriched),
         labels_gp =  gpar(fontsize = 20), border = "black", direction = "horizontal", nrow = 1),
  Legend(col_fun = colors_heatmapbody, 
         title = "snRNA expression", 
         title_gp = gpar(fontsize = 20),
         labels_gp = gpar(fontsize = 20),
         legend_width = unit(6, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_emtscores, 
         title = "EMT score", 
         title_gp = gpar(fontsize = 20),
         labels_gp = gpar(fontsize = 20),
         legend_width = unit(7, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"))

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na, border = "black",
                             ## row
                             show_row_names = T, row_names_gp = gpar(fontsize = 20, fontface = "italic"), row_names_side = "left",
                             row_split = row_split_factor,
                             row_title_rot = 0, row_title_gp = gpar(fontsize = 20),
                             # row_labels = factor_cellgroup,
                             cluster_row_slices = F, show_row_dend = F, 
                             ## column
                             show_column_dend = F, cluster_columns = F, 
                             column_order = column_order_vec,
                             top_annotation = colanno_obj, show_column_names = T, column_title = NA,
                             show_heatmap_legend = F)
file2write <- paste0(dir_out, "EMT_Genes_by_tumorcluster", ".png")
png(file2write, width = 3000, height = 1500, res = 150)
draw(object = p, 
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()

p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na, border = "black",
                             ## row
                             show_row_names = T, row_names_gp = gpar(fontsize = 20, fontface = "italic"), row_names_side = "left",
                             row_split = row_split_factor,
                             row_title_rot = 0, row_title_gp = gpar(fontsize = 20),
                             # row_labels = factor_cellgroup,
                             cluster_row_slices = F, show_row_dend = F, 
                             ## column
                             show_column_dend = F, cluster_columns = F, 
                             column_order = column_order_vec,
                             top_annotation = colanno_obj, show_column_names = F, column_title = NA,
                             show_heatmap_legend = F)
file2write <- paste0(dir_out, "EMT_Genes_by_tumorcluster", ".pdf")
pdf(file2write, width = 20, height = 6.5, useDingbats = F)
draw(object = p,
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()

