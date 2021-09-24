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
avgexp_df <- fread(input = "./Resources/snATAC_Processed_Data/Gene_Activity/AverageGeneActivity_ByTumorPTLOHSubcluster.v2.20210917.tsv", data.table = F)
# ## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/gene_activity/calculate_scores/calculate_msigdb_geneset_scores_wgeneactivity/20210921.v1/MSigDB.Hallmark.tsv")
## input score pre-calculated
epi_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/gene_activity/calculate_scores/calculate_Epithelial_scores_EMTmoduledown_wPT_wgeneactivity/20210921.v1/EpithelialScore.tsv")
s12_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/gene_activity/calculate_scores/calculate_PTS12_scores_wPT_wgeneactivity/20210921.v1/PTS12Score.tsv")
s3_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/gene_activity/calculate_scores/calculate_PTS3_scores_wPT_wgeneactivity/20210921.v1/PTS3Score.tsv")

# preprocess-------------------------------------------------
## specify genes to filter 
emt_genes_df <- data.frame(gene = c("VIM", "FN1", "CDH2", "SERPINE1", "TGFBI",
                                    "LRP2", "CUBN", "SLC17A3", "GATM", "GLYAT",
                                    "SLC5A12", "SLC5A2", 
                                    "SLC13A3", 
                                    "SLC3A1", "SLC16A9", "SLC38A3"),
                           Text_Gene_Group = c(rep("Mesenchymal\nmarkers", 5),
                                               rep("Epithelial/\nproximal-tubule\nmarkers", 5),
                                               rep("PT S1/2", 3), rep("PT S3", 3)))
## add name for the marker groups
genes2filter <- emt_genes_df$gene

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycluster_byaliquot = gsub(x = variable, pattern = "ATACGeneActivity\\.", replacement = "")) %>%
  dplyr::filter(!(grepl(x = id_bycluster_byaliquot, pattern = "LOH"))) %>%
  dplyr::filter(!(grepl(x = id_bycluster_byaliquot, pattern = "NA"))) %>%
  # dplyr::filter((id_bycluster_byaliquot %in% enrich_df$cluster_name) | (grepl(x = id_bycluster_byaliquot, pattern = "PT"))) %>%
  mutate(easyid_column = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(cluster_name = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,2])
## filter out non-tumor and NA tumor cluster
plot_data_long_df <- plot_data_long_df %>%
  filter(!(cluster_name %in% c("", "CNA"))) %>%
  filter(!(grepl(x = id_bycluster_byaliquot, pattern = "00359")))
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
# colors_heatmapbody = colorRamp2(c(-1, 
#                                   0, 
#                                   1), 
#                                 c(color_blue, "white", color_red))
colors_heatmapbody = colorRamp2(c(-1.5, 
                                  0, 
                                  1.5), 
                                c(color_blue, "white", color_red))
## make colors for epithelial score
# summary(scores_epithelial)
color_brown <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[7]
color_green <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[3]
colors_scores_epithelial = colorRamp2(c(100, 0, -100), 
                                      c(brewer.pal(n = 9, name = "YlGn")[c(8,3)], "white"))
## colors for EMT group
colors_emtgroups <- c("black", "white")
names(colors_emtgroups) <- c("High", "Low")
colors_emtenriched <- c("black", "white smoke")
names(colors_emtenriched) <- c("TRUE", "FALSE")
colors_emtscores <- colorRamp2(c(-100, 0, 100), 
                               c("white", brewer.pal(n = 11, name = "Spectral")[c(6,3)]))
# summary(scores_s12)
colors_s12score <- colorRamp2(c(-67.47553, 0, 67.47553), 
                              c("white", brewer.pal(n = 9, name = "BuPu")[c(3,8)]))
# summary(scores_s3)
colors_s3score <- colorRamp2(c(-59.14322, 0, 59.14322), 
                             c("white", brewer.pal(n = 11, name = "BrBG")[c(6,2)]))
## make colors for cell type
colors_celltype <- colors_cellgroup14[c("Tumor cells", "Normal epithelial cells")]
names(colors_celltype) <- c("Tumor cells", "Proximal tubule")
## make colors for S1/2/3 group
colors_s123_group <- c(brewer.pal(n = 9, name = "BuPu")[c(8)], brewer.pal(n = 11, name = "BrBG")[c(2)], 
                       brewer.pal(n = 11, name = "PuOr")[c(11)], brewer.pal(n = 8, name = "Set2")[c(8)])
names(colors_s123_group) <- c("S1/2 enriched", "S3 enriched", "mixed S1/2/3\nidentity", "weak segmental\nidentity")

# get row/column ids ----------------------------------------------------------
columnnames_plot <- colnames(plot_data_mat)
ids_aliquot_wu <- str_split_fixed(string = columnnames_plot, pattern = "_", n = 2)[,1]; ids_aliquot_wu
ids_cluster <- str_split_fixed(string = columnnames_plot, pattern = "_", n = 2)[,2]; ids_cluster
rownames_plot <- rownames(plot_data_mat)

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = rownames_plot, from = emt_genes_df$gene, to = as.vector(emt_genes_df$Text_Gene_Group))
row_split_vec
row_split_factor <- factor(x = row_split_vec, levels = c("Mesenchymal\nmarkers", "Epithelial/\nproximal-tubule\nmarkers", "Tumor-cell\nmarkers", "PT S1/2", "PT S3"))


# make column annotation --------------------------------------------------
## make
emt_scores_vec <- mapvalues(x = columnnames_plot, from = enrich_df$cluster_name, to = as.vector(enrich_df$EPITHELIAL_MESENCHYMAL_TRANSITION_Score)); emt_scores_vec <- as.numeric(emt_scores_vec)
scores_epithelial <- mapvalues(x = columnnames_plot, from = epi_scores_df$cluster_name, to = as.vector(epi_scores_df$score)); scores_epithelial <- as.numeric(scores_epithelial)
scores_s12 <- mapvalues(x = columnnames_plot, from = s12_scores_df$cluster_name, to = as.vector(s12_scores_df$score)); scores_s12 <- as.numeric(scores_s12)
scores_s3 <- mapvalues(x = columnnames_plot, from = s3_scores_df$cluster_name, to = as.vector(s3_scores_df$score)); scores_s3 <- as.numeric(scores_s3)
## make epithelial enriched
emt_enriched_vec <- as.character(columnnames_plot %in% enrich_df$cluster_name[enrich_df$EMT])
## make data frame
colanno_df <- data.frame(columnname = columnnames_plot,
                         cell_type = ifelse(grepl(pattern = "PT", x = columnnames_plot), "Proximal tubule", "Tumor cells"),
                         epithelial_score = scores_epithelial,
                         epithelial_group = ifelse(scores_epithelial >= quantile(x = scores_epithelial, 0.7), "Epi.-strong",
                                                   ifelse(scores_epithelial <= quantile(x = scores_epithelial, 0.3),
                                                          ifelse(columnnames_plot %in% enrich_df$cluster_name[enrich_df$EMT], "EMT", "Epi.-weak"), "Epi.-intermediate")),
                         s123_group = ifelse(scores_s12 >= quantile(scores_s12, 0.9),
                                             ifelse(scores_s3 >=  quantile(scores_s3, 0.9), "mixed S1/2/3\nidentity", "S1/2 enriched"),
                                             ifelse(scores_s3 >=  quantile(scores_s3, 0.9), "S3 enriched", "weak segmental\nidentity")))
## make highlighted samples
index_highlight <- which(columnnames_plot %in% c("C3L.00079.T1_C4", "C3L.00079.T1_C1",  "C3L.00079.T1_C2",  "C3L.00079.T1_C3",  "C3L.00079.T1_C4",
                                                 "C3N.01200.T2_C1", "C3N.01200.T2_C2",
                                                 "C3N.00242.T1_C1"))
texts_highlight <- columnnames_plot[index_highlight];
## make column annotation object
colanno_obj = HeatmapAnnotation(#link = anno_mark(at = index_highlight, labels = texts_highlight, labels_gp = gpar(fontsize = 15), side = "top"),
  CellType = anno_simple(x = colanno_df$cell_type, col = colors_celltype[colanno_df$cell_type], height = unit(0.5, "cm")),
  EpithelialScore = anno_simple(x = colanno_df$epithelial_score, col = colors_scores_epithelial, height = unit(0.75, "cm")),
  # EMTEnriched = anno_simple(x = emt_enriched_vec, col = colors_emtenriched[emt_enriched_vec], height = unit(1, "cm")),
  EMTScore = anno_simple(x = emt_scores_vec, col = colors_emtscores, height = unit(0.75, "cm")),
  PT_S1_S2_Score = anno_simple(x = scores_s12, col = colors_s12score, height = unit(0.75, "cm")),
  PT_S3_Score = anno_simple(x = scores_s3, col = colors_s3score, height = unit(0.75, "cm")),
  PTSegmentSignature = anno_simple(x = colanno_df$s123_group, col = colors_s123_group[colanno_df$s123_group], height = unit(0.75, "cm")),
  annotation_name_gp = gpar(fontsize = 20, fontface = "bold"), annotation_name_side = "left")

# make column order --------------------------------------------------
column_order_vec <- order(scores_epithelial, decreasing = T)

# make column split -------------------------------------------------------
column_group_vec <- colanno_df$epithelial_group
column_group_factor <- factor(x = column_group_vec, levels = c("Epi.-strong", "Epi.-intermediate", "Epi.-weak", "EMT"))

# make legend list --------------------------------------------------------
list_lgd = packLegend(
  # Legend(title = "expression signature enriched", 
  #        title_gp = gpar(fontsize = 20),
  #        labels = names(colors_emtenriched), legend_gp = gpar(fill = colors_emtenriched),
  #        labels_gp =  gpar(fontsize = 20), border = "black", direction = "horizontal", nrow = 1),
  Legend(title = "Cell type", 
         title_gp = gpar(fontsize = 20),
         labels = names(colors_celltype), legend_gp = gpar(fill = colors_celltype), grid_width = unit(0.5, "cm"), grid_height = unit(0.75, "cm"),
         labels_gp =  gpar(fontsize = 20),direction = "horizontal", nrow = 1),
  Legend(col_fun = colors_scores_epithelial, 
         title = "Epithelial score", 
         title_gp = gpar(fontsize = 20),
         labels_gp = gpar(fontsize = 20),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_emtscores, 
         title = "EMT score", 
         title_gp = gpar(fontsize = 20),
         labels_gp = gpar(fontsize = 20),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_s12score, 
         title = "PT S1/S2 score", 
         title_gp = gpar(fontsize = 20),
         labels_gp = gpar(fontsize = 20),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_s3score, 
         title = "PT S3 score", 
         title_gp = gpar(fontsize = 20),
         labels_gp = gpar(fontsize = 20),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),
  Legend(title = "Proximal tubule S1/S2/S3 signature", 
         title_gp = gpar(fontsize = 20),
         labels = names(colors_s123_group), legend_gp = gpar(fill = colors_s123_group), grid_width = unit(0.5, "cm"), grid_height = unit(0.75, "cm"),
         labels_gp =  gpar(fontsize = 20), direction = "horizontal", nrow = 1),
  Legend(col_fun = colors_heatmapbody, 
         title = "Gene activity", 
         title_gp = gpar(fontsize = 20),
         labels_gp = gpar(fontsize = 20),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"),  column_gap = unit(18, "mm"), direction = "horizontal", max_width = unit(36, "cm"))

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
                             show_column_dend = F, cluster_columns = F, cluster_column_slices = F,
                             column_order = column_order_vec, 
                             column_split = column_group_factor, column_title_gp = gpar(fontsize = 25),
                             top_annotation = colanno_obj,
                             show_column_names = T, column_title = NA,
                             show_heatmap_legend = F)
file2write <- paste0(dir_out, "EMT_Genes_by_tumorcluster", ".png")
png(file2write, width = 2500, height = 1500, res = 150)
draw(object = p, 
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()

p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat,
                             col = colors_heatmapbody,
                             na_col = color_na, #border = "black",
                             ## row
                             show_row_names = T, row_names_gp = gpar(fontsize = 20, fontface = "italic"), row_names_side = "left",
                             row_split = row_split_factor,
                             row_title_rot = 0, row_title_gp = gpar(fontsize = 20),
                             # row_labels = factor_cellgroup,
                             cluster_row_slices = F, show_row_dend = F,
                             ## column
                             show_column_dend = F, cluster_columns = F, cluster_column_slices = F,
                             column_order = column_order_vec, column_split = column_group_factor, column_title_gp = gpar(fontsize = 23),
                             top_annotation = colanno_obj,
                             show_column_names = F, column_title = NA,
                             show_heatmap_legend = F)
file2write <- paste0(dir_out, "EMT_Genes_by_tumorcluster", ".pdf")
pdf(file2write, width = 10, height = 8.5, useDingbats = F)
draw(object = p,
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()


