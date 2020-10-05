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
version_tmp <- 2
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
emt_genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/combine_pt_with_emt_markers_all/20200920.v1/Kidney_Specific_EMT_Genes.20200920.v1.tsv")
## input the id meta data
idmetada_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# get genes to plot -------------------------------------------------------
genes_plot <- emt_genes_df$Gene[emt_genes_df$Gene_Group2 == "Mesenchymal"]
genes_plot <- c("FN1", "CDH2", "VIM")

# make ids to plot --------------------------------------------------------
# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
# tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot <- idmetada_df$Aliquot.bulk[idmetada_df$snRNA_available & !is.na(idmetada_df$Aliquot.bulk) & idmetada_df$Sample_Type == "Tumor"]
tumor_bulk_aliquot_ids2plot

case_ids2plot <- mapvalues(x = tumor_bulk_aliquot_ids2plot, from = bulk_meta_tab$Specimen.Label, to = bulk_meta_tab$Case.ID)
case_ids2plot

# normal_bulk_aliquot_ids2plot <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal"]
normal_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Normal" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Normal" & bulk_meta_tab$Set.A == "yes"]))
normal_bulk_aliquot_ids2plot <- normal_bulk_aliquot_ids2plot[normal_bulk_aliquot_ids2plot != case_ids2plot]
normal_bulk_aliquot_ids2plot

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
mat2plot <- mat2plot[genes_plot[genes_plot %in% genes_filtered],]

# make ids ----------------------------------------------------------------
ids_aliquot <- colnames(mat2plot)
sampletypes <- ifelse(ids_aliquot %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal")
labels_col <- mapvalues(x = ids_aliquot, from = bulk_meta_tab$Specimen.Label, to = as.vector(bulk_meta_tab$Case.ID))
labels_col

# make colors -------------------------------------------------------------
colors_sampletype <- c("Tumor" = "red", "Normal" = "green")
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(mat2plot))
heatmapbody_col_fun = colorRamp2(c(-1, 
                                   0, 
                                   1), 
                                 c(color_blue, "white", color_red))

# make column annotation --------------------------------------------------
ca = HeatmapAnnotation(Sample_Type = sampletypes,
                       col = list(Sample_Type = colors_sampletype), 
                       show_legend = F)


# plot heatmap ------------------------------------------------------------
p <- Heatmap(mat2plot, col = heatmapbody_col_fun,
             ## row
             cluster_rows = T, show_row_dend = F,
             row_names_gp = gpar(fontsize = 12), row_names_side = "left",
             ## column
             column_split = sampletypes,
             cluster_columns = T, show_column_dend = F, 
             column_title_gp = gpar(fontsize = 18),
             show_column_names = T, column_labels = labels_col,
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
png(file2write, width = 1000, height = 500, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)
dev.off()


