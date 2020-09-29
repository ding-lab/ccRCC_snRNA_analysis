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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input summary for cell-cell interactin
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/rank_across_celltypes_by_pair_by_avgsigmean/20200925.v3/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/format_expression/combine_avexp_sct_usescale_immunegroups_with_stroma_and_epithelium/20200908.v1/Combined.Average_Expression.SCT.UseScale.20200908.v1.tsv", data.table = F)
## input celltypes-to-cellgroups
celltypes2cellgroup_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/format_expression/combine_avexp_sct_usescale_immunegroups_with_stroma_and_epithelium/20200908.v1/Combined.CellTypes2CellGroup.20200908.v1.tsv")
## input complex annotation
complex2gene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/annotate_interactions/annotate_cellphonedb_complex_to_genesymbols/20200928.v1/CellPhoneDB_complexs_genesymbols.tsv")

# specify pairs to filter -------------------------------------------------
paired_cellgroups.general_process <- "Tumor&Stroma"
summary_filtered_df <- summary_df %>%
  filter(paired_cellgroups.general == paired_cellgroups.general_process) %>%
  filter(rank_byinteraction_acrosscelltypes == 1) %>%
  mutate(paired_celltypes_group = paste0("Tumor&",ifelse(Cell_type.source == "Tumor cells", Cell_type.target, Cell_type.source))) %>%
  arrange(desc(avg_sig_mean), number_sig_cases)
## extract the receptor names
receptors_filter <- unique(summary_filtered_df$gene.target)
receptor_genes <- mapvalues(x = receptors_filter, from = complex2gene_df$identifiers_CellPhoneDB, to = as.vector(complex_cp_df$components_genesymbols))
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

# make split for cell types ------------------------------------------
col_cellgroups_vec <- mapvalues(x = celltypes_plot, from = celltypes2cellgroup_df$colname_celltype, to = as.vector(celltypes2cellgroup_df$Cell_group.shorter))
col_cellgroups_factor <- factor(x = col_cellgroups_vec, levels = c("Nephron_Epithelium", "Stroma", "Immune"))

# make split for genes -----------------------------------------------------
row_genegroup_vec <- as.vector(receptor2gene_df$receptor_name[genes_filter %in% plot_data_df$gene])
row_genegroup_factor <- factor(x = row_genegroup_vec, levels = unique(receptor2gene_df$receptor_name))

# plot heatmap body -----------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = col_fun, na_col = color_na, border = "black",
                             ## row
                             show_row_names = T, row_names_side = "left",
                             show_row_dend = F, cluster_rows = F,
                             row_split = row_genegroup_factor, row_title_rot = 0, cluster_row_slices = F,
                             ## column
                             show_column_names = T, column_names_side = "top",
                             column_names_gp = gpar(fontsize = 15, fontface = "bold"),
                             show_column_dend = F, 
                             column_split = col_cellgroups_factor, column_title = NULL, cluster_column_slices = F,
                             show_heatmap_legend = F)
p
## save heatmap as png
png(filename = paste0(dir_out, "heatmap", ".png"), 
    width = 800, height = 1600, res = 150)
draw(object = p)
dev.off()


