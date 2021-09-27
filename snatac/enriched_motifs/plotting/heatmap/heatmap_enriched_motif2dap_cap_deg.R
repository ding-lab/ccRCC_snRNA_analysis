# Yige Wu @WashU May 2021

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
## input motif mapped
motifs_all_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Peak_Annotation/Motifs_matched.DEG_associated_Peaks.Motif_annotation.20210517.v1.tsv")
## input pathway annotation
# gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/map_genes2pathway__tumorcells_up_degs_msigdb_H_CP_sig_pathways/20210430.v1/DEG2Pathway.tsv")
gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/map_genes2pathway_tumorcells_up_degs_msigdb_H_CP/20210430.v1/DEG2Pathway.tsv")

# preprocess --------------------------------------------------------------
## decide on the motifs to plot
motifs_plot <- c("NFKB2", "NFKB1", "HIF1A", "ARNT::HIF1A", "RBPJ", "MXI1", "ZNF75D", "HSF2", "NEUROD1", "SREBF2", "NEUROG2(var.2)",
                 "KLF15", "NRF1", "SP9", "ZBTB14", "EGR1", "SP3", "TCFL5", "ZNF148", "KLF14", "SP1")
# motifs_plot <- c("NFKB2", "NFKB1", "HIF1A", "ARNT::HIF1A", "RBPJ", "MXI1", "ZNF75D", "HSF2", "NEUROD1", "SREBF2", "NEUROG2(var.2)", "RELA", "MZF1", "MAZ", "SMAD2::SMAD3::4", "ZNF410")
## filter motif-mapping result
motifs_df <- motifs_all_df %>%
  # filter((Peak_Type == "Enhancer") | (Peak_Type == "Promoter")) %>%
  mutate(TF_name = motif.name) %>%
  filter((Peak_Type == "Enhancer") | (Peak_Type == "Promoter" & Motif_Type == "Promoter")) %>%
  filter(TF_name %in% motifs_plot)

## decide on the pathways to plot
pathways_filtered_df <- gene2pathway_df %>%
  dplyr::filter(GeneSymbol %in% motifs_df$Gene[motifs_df$TF_name %in% motifs_plot]) %>%
  select(GeneSet_Name) %>%
  table() %>%
  as.data.frame() %>%
  rename(GeneSet_Name = '.') %>%
  arrange(desc(Freq))
geneset_names_plot <- head(x = pathways_filtered_df$GeneSet_Name, n = 20); geneset_names_plot <- as.vector(geneset_names_plot)
geneset_names_plot <- geneset_names_plot[!(geneset_names_plot %in% c("REACTOME_METABOLISM_OF_CARBOHYDRATES", "WP_PHOTODYNAMIC_THERAPYINDUCED_HIF1_SURVIVAL_SIGNALING", "KEGG_PATHWAYS_IN_CANCER", "REACTOME_SIGNALING_BY_NUCLEAR_RECEPTORS", "WP_NUCLEAR_RECEPTORS_METAPATHWAY"))]
# geneset_names_plot <- c(geneset_names_plot[geneset_names_plot != "WP_PATHWAYS_IN_CLEAR_CELL_RENAL_CELL_CARCINOMA"], "WP_PATHWAYS_IN_CLEAR_CELL_RENAL_CELL_CARCINOMA")
geneset_names_plot <- c(geneset_names_plot[geneset_names_plot != "WP_PATHWAYS_IN_CLEAR_CELL_RENAL_CELL_CARCINOMA"], "WP_PATHWAYS_IN_CLEAR_CELL_RENAL_CELL_CARCINOMA", "REACTOME_SIGNALING_BY_NOTCH", "HALLMARK_NOTCH_SIGNALING")
## decide on the genes to plot
genes_pathwayselect <- gene2pathway_df$GeneSymbol[gene2pathway_df$GeneSet_Name %in% geneset_names_plot]; genes_pathwayselect <- unique(genes_pathwayselect)
genes_plot_df <- motifs_df %>%
  filter(TF_name %in% motifs_plot) %>%
  filter(Gene %in% genes_pathwayselect) %>%
  select(Gene) %>%
  unique()
# genes_plot <- head(x = genes_plot_df$Gene, n = 50)
genes_plot <- head(x = genes_plot_df$Gene, n = 120)
## make column order
gene2pathway_filtered_df <- gene2pathway_df %>%
  filter(GeneSymbol %in% genes_plot) %>%
  filter(GeneSet_Name %in% geneset_names_plot)
gene2pathway_filtered_df$GeneSet_Order <- mapvalues(x = gene2pathway_filtered_df$GeneSet_Name, from = pathways_filtered_df$GeneSet_Name, to = as.vector(pathways_filtered_df$Freq))
gene2pathway_filtered_df$GeneSet_Order <- as.numeric(gene2pathway_filtered_df$GeneSet_Order)
gene2pathway_filtered_df <- gene2pathway_filtered_df %>%
  arrange(desc(GeneSet_Order), GeneSet_Name)
genes_plot_ordered <- unique(gene2pathway_filtered_df$GeneSymbol)

# make plot data ----------------------------------------------------------
plotdata_df <- motifs_df %>%
  filter(TF_name %in% motifs_plot) %>%
  filter(Gene %in% genes_plot) %>%
  select(TF_name, Gene, Peak_Type) %>%
  unique() %>%
  arrange(TF_name, Gene, Peak_Type) %>%
  group_by(TF_name, Gene) %>%
  summarise(DAP_Type.sum = paste0(Peak_Type, collapse = "&"))
plotdata_wide_df <- dcast(data = plotdata_df, formula = TF_name ~ Gene, value.var = "DAP_Type.sum")
plotdata_mat <- as.matrix(x = plotdata_wide_df[,-1])
rownames(plotdata_mat) <- plotdata_wide_df$TF_name
plotdata_mat <- plotdata_mat[,genes_plot_ordered]
## change values
# plotdata_mat[plotdata_mat == "Enhancer"] <- 1
# plotdata_mat[plotdata_mat == "Enhancer&Promoter"] <- 2
# plotdata_mat[plotdata_mat == "Promoter"] <- 3
# plotdata_mat[is.na(plotdata_mat)] <- 0

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "white"
## specify colors for the heatmap body
colors_heatmapbody <- RColorBrewer::brewer.pal(name = "Dark2", n = 3)
names(colors_heatmapbody) <- c("Promoter", "Enhancer&Promoter", "Enhancer")
# colors for in pathway
colors_isgeneset <- c("black", "grey90")
names(colors_isgeneset) <- c("TRUE", "FALSE")

# make column annotation --------------------------------------------------
colanno_raw_df <- dcast(data = gene2pathway_filtered_df, formula = GeneSymbol ~ GeneSet_Name)
colanno_df <- colanno_raw_df[,-1]
colanno_df[!is.na(colanno_df)] <- "TRUE"
colanno_df[is.na(colanno_df)] <- "FALSE"
rownames(colanno_df) <- colanno_raw_df$GeneSymbol
colanno_df <- colanno_df[genes_plot_ordered, geneset_names_plot]
geneset_names_sim_plot <- tolower(x = str_split_fixed(string = geneset_names_plot, pattern = "_", n = 2)[,2])
colnames(colanno_df) <-geneset_names_sim_plot
## make colors
colors_colanno_list <- list()
for (colname_tmp in colnames(colanno_df)) {
  colors_colanno_list[[colname_tmp]] <- colors_isgeneset
}
colanno_obj = HeatmapAnnotation(df = colanno_df, col = colors_colanno_list,
                                annotation_name_gp = gpar(fontsize = 12), annotation_name_side = "left", show_legend = F)

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plotdata_mat, na_col = color_na,
                             col = colors_heatmapbody,
                             ## row
                             cluster_rows = T,
                             show_row_names = T, row_names_gp = gpar(fontsize = 14), row_names_side = "left",
                             # show_row_dend = F, cluster_row_slices = T, row_labels = rowlabels_plot,
                             ## column
                             show_column_dend = F, cluster_columns = F,
                             top_annotation = colanno_obj,
                             # show_column_names = F, column_names_side = "top", column_names_gp = gpar(fontsize = 5),
                             show_heatmap_legend = F)
p

# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(labels = names(colors_heatmapbody), labels_gp = gpar(fontsize = 14),
         title = "Motif location to gene", title_gp = gpar(fontsize = 14),
         legend_gp = gpar(fill = colors_heatmapbody), border = NA),
  Legend(labels = names(colors_isgeneset), labels_gp = gpar(fontsize = 14),
         title = "Gene in pathway", title_gp = gpar(fontsize = 14),
         legend_gp = gpar(fill = colors_isgeneset), border = NA))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Motif2DEG", ".png")
# png(file2write, width = 2000, height = 1000, res = 150)
# png(file2write, width = 3000, height = 1000, res = 150)
png(file2write, width = 3000, height = 1000, res = 150)
draw(object = p, 
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()
file2write <- paste0(dir_out, "Motif2DEG", ".pdf")
# pdf(file2write, width = 20, height = 6.5, useDingbats = F)
pdf(file2write, width = 25, height = 12, useDingbats = F)
draw(object = p,
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()
file2write <- paste0(dir_out, "DEGs.SelectedPathways.tsv")
write.table(x = genes_plot_df, file = file2write, quote = F, sep = "\t", row.names = F)

