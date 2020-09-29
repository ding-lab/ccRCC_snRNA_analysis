# Yige Wu @WashU Jul 2020

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
## input protein data
protein_tab <- fread("~/Box/Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/proteome/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)
## input bulk meta data for the entire set
bulk_meta_tab <- fread("~/Box/Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
## input the peaks to TF motifs
deg2tfgene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_tumor_vs_normal_degs_to_tf_by_snatac/20200916.v1/DEGs_with_DA_peaks.20200916.v1.tsv")

# get genes to plot -------------------------------------------------------
gene_tf <- "NR3C1"
genes_targets <- deg2tfgene_df$target_genesymbol[deg2tfgene_df$source_genesymbol == gene_tf]
genes_targets <- genes_targets[genes_targets != "GATM"]
genes_plot <- c(gene_tf, genes_targets)

# make ids to plot --------------------------------------------------------
# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
normal_bulk_aliquot_ids2plot <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal"]
normal_bulk_aliquot_ids2plot

case_ids2plot <- mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$Specimen.Label, to = bulk_meta_tab$Case.ID)
case_ids2plot

tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot

# make heatmap body -------------------------------------------------------
## make the matrix to plot the heatmap
protein_tab2plot <- protein_tab %>%
  filter(Index %in% genes_plot) %>%
  select("Index", "ReferenceIntensity", tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)

protein_mat2plot <- protein_tab2plot %>%
  select(-Index) %>%
  select(-ReferenceIntensity)
rownames(protein_mat2plot) <- protein_tab2plot$Index

mat2plot <- as.matrix(protein_mat2plot) - as.vector(protein_tab2plot$ReferenceIntensity)
genes_filtered <- rowSums(!is.na(mat2plot))
genes_filtered <- names(genes_filtered)[genes_filtered > 0]
genes_filtered <- c("NR3C1", genes_filtered[!(genes_filtered %in% c("NR3C1", "GATM"))])
genes_filtered <- c("NR3C1", "PCSK6", "KCTD3", "EGFR", "PDIA5", "SEMA6A", "PDK4", "COL23A1", "SEMA5B", "PLOD2")
mat2plot <- mat2plot[genes_filtered,]
ids_aliquot <- colnames(mat2plot)
sampletypes <- ifelse(ids_aliquot %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal")
# make colors -------------------------------------------------------------
colors_sampletype <- c("Tumor" = "red", "Normal" = "green")
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
heatmapbody_col_fun = colorRamp2(c(-2, 
                                   0, 
                                   2), 
                                 c(color_blue, "white", color_red))

# make column annotation --------------------------------------------------
ca = HeatmapAnnotation(Sample_Type = sampletypes,
                       col = list(Sample_Type = colors_sampletype), 
                       show_legend = F)


# plot heatmap ------------------------------------------------------------
p <- Heatmap(mat2plot, col = heatmapbody_col_fun,
             # top_annotation = ca,
             column_split = sampletypes,
             cluster_columns = T, show_column_dend = F, 
             column_title_gp = gpar(fontsize = 18),
             cluster_rows = F, show_row_dend = F,
             row_names_gp = gpar(fontsize = 12), row_names_side = "left",
             show_column_names = F, 
             show_heatmap_legend = F)

annotation_lgd = list(
  Legend(col_fun = heatmapbody_col_fun, 
         title = "Protein abundance value\nlog2Intensity\n(Sample-Reference)", 
         legend_width = unit(6, "cm"),
         direction = "horizontal"))


## save heatmap to file
file2write <- paste0(dir_out, "Heatmap_by_sample.", ".pdf")
pdf(file2write, width = 20, height = 5)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)
dev.off()
file2write <- paste0(dir_out, "Heatmap_by_sample.",".png")
png(file2write, width = 1200, height = 500, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)
dev.off()


