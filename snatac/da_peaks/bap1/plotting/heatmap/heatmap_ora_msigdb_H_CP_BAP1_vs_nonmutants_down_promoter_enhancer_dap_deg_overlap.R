# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
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

# input data --------------------------------------------------------------
ora_result_df1 <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/pathway/ora_msigdb_H_CP_BAP1_vs_nonmutants_down_promoter_enhancer_dap_deg_overlap/20210913.v2/ORA_Results.tsv")
ora_result_df2 <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/pathway/ora_msigdb_H_CP_BAP1_vs_nonmutants_up_promoter_enhancer_dap_deg_overlap/20210913.v1/ORA_Results.tsv")
dap2deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/overlap_degs/overlap_bap1_vs_nonmutant_enhancer_promoter_peaks_with_degs/20210913.v1/BAP1_vs_NonMutant_DAP2DEG.20210913.v1.tsv")

# select gene sets --------------------------------------------------------
enricher_filtered_df <- ora_result_df1 %>%
  filter(pvalue < 0.05) %>%
  filter(Count >= 2) %>%
  mutate(row_id = ID)
enricher_filtered_df$Keep <- F
## step one: add one pathway at a time based on the max overlap
### for each group, add the top pathway first and add one for each iteration with % overlap < 50%
enricher_filtered_df$Keep[enricher_filtered_df$ID %in% c("KEGG_ADHERENS_JUNCTION", "REACTOME_RHOA_GTPASE_CYCLE", "KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM", "WP_WNT_SIGNALING",
                                                         "HALLMARK_KRAS_SIGNALING_DN", "HALLMARK_P53_PATHWAY", "REACTOME_BIOLOGICAL_OXIDATIONS", "HALLMARK_ESTROGEN_RESPONSE_EARLY",
                                                         "REACTOME_TNFS_BIND_THEIR_PHYSIOLOGICAL_RECEPTORS", "WP_AMINO_ACID_METABOLISM", 
                                                         "REACTOME_CARGO_CONCENTRATION_IN_THE_ER", "WP_NUCLEAR_RECEPTORS_METAPATHWAY")] <- T

## step 2: run this
enricher_filtered_df$max_overlap_ratio <- sapply(enricher_filtered_df$row_id, function(row_id_tmp, test_df) {
  # sapply(head(enricher_filtered_df$row_id), function(row_id_tmp, test_df) {
  
  test_tmp <- test_df$test[test_df$row_id == row_id_tmp]
  pvalue_tmp <- test_df$pvalue[test_df$row_id == row_id_tmp]
  id_tmp <- test_df$ID[test_df$row_id == row_id_tmp]
  test_keep_df <- test_df[test_df$Keep & test_df$pvalue <= pvalue_tmp & test_df$ID != id_tmp,]
  
  genes_tmp <- str_split(string = test_df$geneID[test_df$row_id == row_id_tmp], pattern = "\\/")[[1]]
  max_overlap_ratio <- 0
  for (genestring_kept_tmp in test_keep_df$geneID) {
    genes_kept_tmp <- str_split(string = genestring_kept_tmp, pattern = "\\/")[[1]]
    genes_common_tmp <- intersect(genes_tmp, genes_kept_tmp)
    overlap_ratio_tmp <- length(genes_common_tmp)/length(genes_tmp)
    max_overlap_ratio <- max(c(max_overlap_ratio, overlap_ratio_tmp))
  }
  return(max_overlap_ratio)
}, test_df = enricher_filtered_df)
## go back to step 1


enricher_filtered_df <- ora_result_df2 %>%
  filter(pvalue < 0.05) %>%
  filter(Count >= 2) %>%
  mutate(row_id = ID)
enricher_filtered_df$Keep <- F
## step one: add one pathway at a time based on the max overlap
### for each group, add the top pathway first and add one for each iteration with % overlap < 50%
enricher_filtered_df$Keep[enricher_filtered_df$ID %in% c("WP_PATHWAYS_IN_CLEAR_CELL_RENAL_CELL_CARCINOMA", "REACTOME_EPHA_MEDIATED_GROWTH_CONE_COLLAPSE",
                                                         "WP_BRAINDERIVED_NEUROTROPHIC_FACTOR_BDNF_SIGNALING_PATHWAY", "HALLMARK_HEME_METABOLISM")] <- T

## step 2: run this
enricher_filtered_df$max_overlap_ratio <- sapply(enricher_filtered_df$row_id, function(row_id_tmp, test_df) {
  # sapply(head(enricher_filtered_df$row_id), function(row_id_tmp, test_df) {
  
  test_tmp <- test_df$test[test_df$row_id == row_id_tmp]
  pvalue_tmp <- test_df$pvalue[test_df$row_id == row_id_tmp]
  id_tmp <- test_df$ID[test_df$row_id == row_id_tmp]
  test_keep_df <- test_df[test_df$Keep & test_df$pvalue <= pvalue_tmp & test_df$ID != id_tmp,]
  
  genes_tmp <- str_split(string = test_df$geneID[test_df$row_id == row_id_tmp], pattern = "\\/")[[1]]
  max_overlap_ratio <- 0
  for (genestring_kept_tmp in test_keep_df$geneID) {
    genes_kept_tmp <- str_split(string = genestring_kept_tmp, pattern = "\\/")[[1]]
    genes_common_tmp <- intersect(genes_tmp, genes_kept_tmp)
    overlap_ratio_tmp <- length(genes_common_tmp)/length(genes_tmp)
    max_overlap_ratio <- max(c(max_overlap_ratio, overlap_ratio_tmp))
  }
  return(max_overlap_ratio)
}, test_df = enricher_filtered_df)
## go back to step 1


# set plotting parameters -------------------------------------------------
genesets_plot <- c("KEGG_ADHERENS_JUNCTION", "REACTOME_RHOA_GTPASE_CYCLE", "KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM", "WP_WNT_SIGNALING",
                   "HALLMARK_KRAS_SIGNALING_DN", "HALLMARK_P53_PATHWAY", "REACTOME_BIOLOGICAL_OXIDATIONS", "HALLMARK_ESTROGEN_RESPONSE_EARLY",
                   "REACTOME_TNFS_BIND_THEIR_PHYSIOLOGICAL_RECEPTORS", "WP_AMINO_ACID_METABOLISM", "WP_NUCLEAR_RECEPTORS_METAPATHWAY",
                   "WP_PATHWAYS_IN_CLEAR_CELL_RENAL_CELL_CARCINOMA", "REACTOME_EPHA_MEDIATED_GROWTH_CONE_COLLAPSE",
                   "WP_BRAINDERIVED_NEUROTROPHIC_FACTOR_BDNF_SIGNALING_PATHWAY", "HALLMARK_HEME_METABOLISM")

# process -----------------------------------------------------------------
ora_result_df <- rbind(ora_result_df1, ora_result_df2)
ora_filtered_df <- ora_result_df %>%
  filter(ID %in% genesets_plot)

genes_list <- sapply(ora_filtered_df$geneID, function(genes_string) {
  genes_vec <- str_split(string = genes_string, pattern = "\\/")[[1]]
  return(genes_vec)
})
genenumber_vec <- sapply(ora_filtered_df$geneID, function(genes_string) {
  genes_vec <- str_split(string = genes_string, pattern = "\\/")[[1]]
  return(length(genes_vec))
})
names(genes_list) <- ora_filtered_df$ID
gene2geneset_long_df <- data.frame(genesymbol = unlist(genes_list), geneset = rep(ora_filtered_df$ID, genenumber_vec))
gene2geneset_long_df$avg_log2FC.snRNA <- mapvalues(x = gene2geneset_long_df$genesymbol, from = dap2deg_df$Gene, to = as.vector(dap2deg_df$avg_log2FC.snRNA))
gene2geneset_long_df$avg_log2FC.snRNA <- as.numeric(gene2geneset_long_df$avg_log2FC.snRNA)

# make heatmap data -------------------------------------------------------
gene2geneset_wide_df <- dcast(data = gene2geneset_long_df, geneset ~ genesymbol, value.var = "avg_log2FC.snRNA")
plotdata_mat <- as.matrix(gene2geneset_wide_df[,-1])
rownames(plotdata_mat) <- gene2geneset_wide_df$geneset
plotdata_mat[is.na(plotdata_mat)] <- 0

# make colors -------------------------------------------------------------
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plotdata_mat))
colors_heatmapbody = colorRamp2(c(-1, 0, 1), c("purple", "grey20", "yellow"))
colors_heatmapbody = colorRamp2(c(-1, 0, 1), c(color_blue, "white", color_red))

fontsize_plot <- 18
list_lgd = packLegend(
  Legend(col_fun = colors_heatmapbody, 
         title = "log2 fold change\nBAP1-mutants vs. non-mutants", 
         title_gp = gpar(fontsize = fontsize_plot),
         labels_gp = gpar(fontsize = fontsize_plot),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"))

# plot --------------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plotdata_mat, col = colors_heatmapbody,
                             cluster_columns = T, show_column_dend = F, column_names_gp = gpar(fontsize = fontsize_plot, fontface = "italic"),
                             cluster_rows = T, show_row_dend = F, row_names_side = "left", row_names_gp = gpar(fontsize = fontsize_plot),
                             show_heatmap_legend = F)
file2write <- paste0(dir_out, "heatmap", ".pdf")
pdf(file2write, width = 14, height = 6.6, useDingbats = F)
draw(object = p,
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

roworder_plot <- rownames(plotdata_mat)[rev(row_order(p))]
columnorder_plot <- colnames(plotdata_mat)[rev(column_order(p))]
plotdata_mat2 <- plotdata_mat[roworder_plot, columnorder_plot]
rowlabels_plot <- str_split_fixed(string = roworder_plot, pattern = "_", n = 2)[,2]
rowlabels_plot[rowlabels_plot == "BRAINDERIVED_NEUROTROPHIC_FACTOR_BDNF_SIGNALING_PATHWAY"] <- "BDNF_SIGNALING_PATHWAY"
rowlabels_plot[rowlabels_plot == "PATHWAYS_IN_CLEAR_CELL_RENAL_CELL_CARCINOMA"] <- "PATHWAY in ccRCC"
p2 <- ComplexHeatmap::Heatmap(matrix = plotdata_mat2, col = colors_heatmapbody,
                             cluster_columns = F,  column_names_gp = gpar(fontsize = fontsize_plot, fontface = "italic"),
                             cluster_rows = F, row_names_side = "left", row_names_gp = gpar(fontsize = fontsize_plot), row_labels = rowlabels_plot,
                             show_heatmap_legend = F)
file2write <- paste0(dir_out, "heatmap2", ".pdf")
pdf(file2write, width = 13, height = 6.6, useDingbats = F)
draw(object = p2,
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()



