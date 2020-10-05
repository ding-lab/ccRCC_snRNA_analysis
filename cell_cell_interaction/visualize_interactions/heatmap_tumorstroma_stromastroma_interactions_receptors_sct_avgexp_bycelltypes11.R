# Yige Wu @WashU Sep 2020

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
## input summary for cell-cell interactin
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/annotate_interactions/annotate_singlegeneexp_annotated_interactions_with_druggability/20200929.v1/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/format_expression/combine_avexp_sct_usescale_immunegroups_with_stroma_and_epithelium/20200908.v1/Combined.Average_Expression.SCT.UseScale.20200908.v1.tsv", data.table = F)
## input celltypes-to-cellgroups
celltypes2cellgroup_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/format_expression/combine_avexp_sct_usescale_immunegroups_with_stroma_and_epithelium/20200908.v1/Combined.CellTypes2CellGroup.20200908.v1.tsv")
## input complex annotation
complex2gene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/annotate_interactions/annotate_cellphonedb_complex_to_genesymbols/20200928.v1/CellPhoneDB_complexs_genesymbols.tsv")
## input bulk average protein expression
avg_bulkpro_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/Tumor_normal_bulk_Protein.txt")
## input bulk average rna expression
avg_bulkrna_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/Tumor_normal_bulk_RNA.txt")

# specify pairs to filter -------------------------------------------------
paired_cellgroups.general_process <- "Tumor&Stroma"
summary_filtered_df1 <- summary_df %>%
  filter(paired_cellgroups.general == paired_cellgroups.general_process) %>%
  filter(rank_byinteraction_acrosscelltypes == 1) %>%
  filter(rank_genesource_acrosscelltypes == 1 & rank_genetarget_acrosscelltypes == 1) %>%
  filter(!(interacting_pair %in% c("TNC_aVb6 complex", "PVR_NECTIN3", "FGF7_FGFR4"))) %>%
  mutate(paired_celltypes_group = paste0("Tumor&",ifelse(Cell_type.source == "Tumor cells", Cell_type.target, Cell_type.source))) %>%
  arrange(desc(avg_sig_mean), number_sig_cases)
paired_cellgroups.general_process <- "Stroma&Stroma"
summary_filtered_df2 <- summary_df %>%
  filter(paired_cellgroups.general == paired_cellgroups.general_process) %>%
  filter(rank_byinteraction_acrosscelltypes == 1) %>%
  filter(rank_genesource_acrosscelltypes == 1 & rank_genetarget_acrosscelltypes == 1) %>%
  filter(!is.na(gene.source.druggable) | !is.na(gene.target.druggable)) %>%
  filter(!is_integrin) %>%
  filter(!(interacting_pair %in% c("PlexinA2_complex1_SEMA3A", "PlexinA1_complex1_SEMA6D", "ACVR_1B2A receptor_INHBA", "EPOR_KITLG", "FGFR1_FGF7", "FGF1_FGFR1"))) %>%
  mutate(paired_celltypes_group = ifelse(Cell_type.target == "Endothelial cells", paste0("Endothelial cells", "&", Cell_type.source), paste0(Cell_type.source, "&", Cell_type.target))) %>%
  arrange(desc(avg_sig_mean), number_sig_cases)
## extract the ligand genes
receptors_filter <- unique(c(summary_filtered_df1$gene.target, summary_filtered_df2$gene.target))
receptors_filter
receptor_genes <- mapvalues(x = receptors_filter, from = complex2gene_df$identifiers_CellPhoneDB, to = as.vector(complex2gene_df$components_genesymbols))
receptor_genes
## split
genes_filter_list <- sapply(X = receptor_genes, FUN = function(x) {
  text_split <- str_split(string = x, pattern = "_")[[1]]
  return(text_split)
})
genes_filter_list
names(genes_filter_list) <- receptors_filter
## get the repeat index
idx_receptor_rep <- sapply(X = receptor_genes, FUN = function(x) {
  text_split <- str_split(string = x, pattern = "_")[[1]]
  return(length(text_split))
})
## make mapping tale
receptor2gene_df <- data.frame(receptor_name = rep(x = receptors_filter, idx_receptor_rep), genesymbol = unlist(genes_filter_list, use.names = T))
## geth ordered gene list to filter
genes_filter <- as.vector(receptor2gene_df$genesymbol)

# format the column names to only aliquot id ------------------------------
## filtr the rows
plot_data_df <- avgexp_df %>%
  filter(gene %in% genes_filter)
## filter the columns and make data matrix
plot_data_mat <- as.matrix(plot_data_df[,-1])
## add gene name
rownames(plot_data_mat) <- plot_data_df$gene
## order by genes
plot_data_mat <- plot_data_mat[genes_filter[genes_filter %in% plot_data_df$gene],]
plot_data_mat %>% head()

# get dimension names -----------------------------------------------------
genes_plot <- rownames(plot_data_mat)
genes_plot
celltypes_plot <- colnames(plot_data_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for baseline expression
col_baselineexp <- colorRamp2(c(0, 2), c("white", "orange"))
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
col_fun = colorRamp2(c(-1.5, 
                       0, 
                       1.5), 
                     c(color_blue, "white", color_red))
## make colors for bulk protein values
summary(avg_bulkpro_df$Mean_tumor)
color_orange <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[5]
color_purple <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[4]
colors_bulkpro <- colorRamp2(c(-1, 
                               0, 
                               1), 
                             c(color_purple, "white", color_orange))
## make colors for bulk rna values
summary(avg_bulkrna_df$Mean_tumor)
colors_bulkrna <- colorRamp2(c(-0.75, 
                               0, 
                               0.75), 
                             c(color_purple, "white", color_orange))


# make split for cell types ------------------------------------------
col_cellgroups_vec <- mapvalues(x = celltypes_plot, from = celltypes2cellgroup_df$colname_celltype, to = as.vector(celltypes2cellgroup_df$Cell_group.shorter))
col_cellgroups_factor <- factor(x = col_cellgroups_vec, levels = c("Nephron_Epithelium", "Stroma", "Immune"))

# make split for genes -----------------------------------------------------
row_genegroup_vec <- as.vector(receptor2gene_df$receptor_name[genes_filter %in% plot_data_df$gene])
row_genegroup_factor <- factor(x = row_genegroup_vec, levels = unique(receptor2gene_df$receptor_name))

# make row annotation -----------------------------------------------------
## make vector for averaged bulk protein
rownames(avg_bulkpro_df) <- avg_bulkpro_df$Gene
mean_protein_tumor <- avg_bulkpro_df[genes_plot, "Mean_tumor"]; mean_protein_tumor
mean_protein_normal <- avg_bulkpro_df[genes_plot, "Mean_Normal"]; mean_protein_normal
## make vector for averaged bulk rna
rownames(avg_bulkrna_df) <- avg_bulkrna_df$Gene
mean_rna_tumor <- avg_bulkrna_df[genes_plot, "Mean_tumor"]; mean_rna_tumor
mean_rna_normal <- avg_bulkrna_df[genes_plot, "Mean_Normal"]; mean_rna_normal
## make row annotation object
rowanno_obj <- rowAnnotation(Avg_Bulk_mRNA = anno_simple(x = cbind(mean_rna_tumor,mean_rna_normal), col = colors_bulkrna, width = unit(x = 1, "cm")),
                             Avg_Bulk_Protein = anno_simple(x = cbind(mean_protein_tumor, mean_protein_normal), col = colors_bulkpro, width = unit(x = 1, "cm")),
                             annotation_name_side = "top", border = T)

# plot heatmap body -----------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = col_fun, na_col = color_na, border = "black",
                             ## row
                             show_row_names = T, row_names_side = "left", row_names_gp = gpar(fontface = "italic", fontsize = 15),
                             show_row_dend = F, cluster_rows = F,
                             row_split = row_genegroup_factor, row_title_rot = 0, cluster_row_slices = F,
                             right_annotation = rowanno_obj,
                             ## column
                             show_column_names = T, column_names_side = "top",
                             column_names_gp = gpar(fontsize = 15),
                             show_column_dend = F, 
                             column_split = col_cellgroups_factor, cluster_column_slices = F, cluster_columns = F,
                             column_title = NULL,
                             show_heatmap_legend = F)
p
## make legend
list_lgd = list(
  Legend(title = "Scaled snRNA expression", title_gp = gpar(fontsize = 12),
         col_fun = colors_heatmapbody, 
         legend_width = unit(3, "cm"),
         direction = "horizontal"),
  Legend(title = "Bulk mRNA z-score", title_gp = gpar(fontsize = 12),
         col_fun = colors_bulkrna, 
         legend_width = unit(3, "cm"),
         direction = "horizontal"),
  Legend(title = "Bulk protein z-score", title_gp = gpar(fontsize = 12),
         col_fun = colors_bulkpro, 
         legend_width = unit(3, "cm"),
         direction = "horizontal"))

## save heatmap as png
png(filename = paste0(dir_out, "heatmap", ".png"), 
    width = 900, height = 1000, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
file2write <- paste0(dir_out, "heatmap", ".pdf")
pdf(file2write, width = 6, height = 7)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
## save with no legend
file2write <- paste0(dir_out, "heatmap.nolegend", ".pdf")
pdf(file2write, width = 6, height = 6)
draw(object = p)
dev.off()


