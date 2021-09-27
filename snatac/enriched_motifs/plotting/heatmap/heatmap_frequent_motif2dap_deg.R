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
version_tmp <- "top.enrichedmotifs"
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input motif mapped
motifs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/summarize_tf_motifs/count_motif_promoter_enhancer_occurance_in_up_dap_deg/20210512.v1/Motif_in_ProEnh_DAP_DEGs.20210512.v1.tsv")
## input pathway annotation
gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/map_genes2pathway__tumorcells_up_degs_msigdb_H_CP_sig_pathways/20210430.v1/DEG2Pathway.tsv")
## input motif count
count_motifs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/summarize_tf_motifs/count_motif_promoter_enhancer_occurance_in_up_dap_deg/20210512.v1/Count_Motif_in_DAP_DEGs.20210512.v1.tsv")
count_motifs_bypathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/summarize_tf_motifs/count_motif_promoter_enhancer_occurance_in_up_dap_deg/20210512.v1/Count_Motif_in_DAP_DEGs.ByPathway.20210512.v1.tsv")

# preprocess --------------------------------------------------------------
## decide on the motifs to plot
motifs_plot <- head(x = count_motifs_df$TF_name, n = 10)
motifs_plot <- head(x = count_motifs_df$TF_name, n = 20)
motifs_plot <- unique(c(motifs_plot, "HIF1A"))
motifs_plot <- c("NFKB2", "NFKB1", "HIF1A", "ARNT::HIF1A", "RBPJ", "MXI1", "ZNF75D", "HSF2", "NEUROD1", "SREBF2", "NEIROG2(var.2)")

## decide on the pathways to plot
pathways_filtered_df <- count_motifs_bypathway_df %>%
  dplyr::filter(TF_name %in% motifs_plot) %>%
  select(GeneSet_Name, Freq) %>%
  unique() %>%
  dplyr::group_by(GeneSet_Name) %>%
  dplyr::slice_max(order_by = Freq, n = 1) %>%
  arrange(desc(Freq))
geneset_names_plot <- c(head(x = pathways_filtered_df$GeneSet_Name, n = 10), "WP_PATHWAYS_IN_CLEAR_CELL_RENAL_CELL_CARCINOMA", "KEGG_RENAL_CELL_CARCINOMA")
geneset_names_plot <- head(x = pathways_filtered_df$GeneSet_Name, n = 20)

## decide on the genes to plot
genes_pathwayselect <- gene2pathway_df$GeneSymbol[gene2pathway_df$GeneSet_Name %in% geneset_names_plot]; genes_plot <- unique(genes_plot)
# genes_plot <- genes_plot[genes_plot %in% motifs_df$genesymbol_deg[motifs_df$TF_name %in% motifs_plot]]
genes_plot_df <- motifs_df %>%
  filter(TF_name %in% motifs_plot) %>%
  filter(genesymbol_deg %in% genes_pathwayselect) %>%
  select(genesymbol_deg, mean_avg_logFC) %>%
  unique() %>%
  arrange(desc(mean_avg_logFC))
# genes_plot <- head(x = genes_plot_df$genesymbol_deg, n = 50)
genes_plot <- head(x = genes_plot_df$genesymbol_deg, n = 121)
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
  filter(genesymbol_deg %in% genes_plot) %>%
  select(TF_name, genesymbol_deg, DAP_Type.strict) %>%
  unique() %>%
  arrange(TF_name, genesymbol_deg, DAP_Type.strict) %>%
  group_by(TF_name, genesymbol_deg) %>%
  summarise(DAP_Type.sum = paste0(DAP_Type.strict, collapse = "&"))
plotdata_wide_df <- dcast(data = plotdata_df, formula = TF_name ~ genesymbol_deg, value.var = "DAP_Type.sum")
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
png(file2write, width = 3000, height = 1200, res = 150)
draw(object = p, 
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()
file2write <- paste0(dir_out, "Motif2DEG", ".pdf")
# pdf(file2write, width = 20, height = 6.5, useDingbats = F)
pdf(file2write, width = 20, height = 8, useDingbats = F)
draw(object = p,
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()
