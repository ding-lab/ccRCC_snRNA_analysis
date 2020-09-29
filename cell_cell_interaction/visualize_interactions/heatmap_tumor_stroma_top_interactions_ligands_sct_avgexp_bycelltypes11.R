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
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/rank_across_celltypes_by_pair_by_avgsigmean/20200925.v3/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/format_expression/combine_avexp_sct_usescale_immunegroups_with_stroma_and_epithelium/20200908.v1/Combined.Average_Expression.SCT.UseScale.20200908.v1.tsv", data.table = F)
## input celltypes-to-cellgroups
celltypes2cellgroup_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/format_expression/combine_avexp_sct_usescale_immunegroups_with_stroma_and_epithelium/20200908.v1/Combined.CellTypes2CellGroup.20200908.v1.tsv")
## input the DEG  list
deg_df <- fread(input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/annotate_markers/annotate_cellgroup4_degs_with_da_promoter/20200916.v1/CellTypeDEGTop50.Annotation.20200916.v1.tsv", data.table = F)

# specify pairs to filter -------------------------------------------------
paired_cellgroups.general_process <- "Tumor&Stroma"
summary_filtered_df <- summary_df %>%
  filter(paired_cellgroups.general == paired_cellgroups.general_process) %>%
  filter(rank_byinteraction_acrosscelltypes == 1) %>%
  mutate(paired_celltypes_group = paste0("Tumor&",ifelse(Cell_type.source == "Tumor cells", Cell_type.target, Cell_type.source))) %>%
  arrange(desc(avg_sig_mean), number_sig_cases)
## extract the ligand genes
genes_filter <- summary_filtered_df$gene.source
# genes_filter[!(genes_filter %in% deg_df$gene)]

# format the column names to only aliquot id ------------------------------
## filtr the rows
plot_data_df <- avgexp_df %>%
  filter(gene %in% genes_filter)
## filter the columns and make data matrix
plot_data_mat <- as.matrix(plot_data_df[,-1])
## add gene name
rownames(plot_data_mat) <- plot_data_df$gene
## order by genes
plot_data_mat <- plot_data_mat[genes_filter,]
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

# make row split for cell types ------------------------------------------
col_cellgroups_vec <- mapvalues(x = celltypes_plot, from = celltypes2cellgroup_df$colname_celltype, to = as.vector(celltypes2cellgroup_df$Cell_group.shorter))
col_cellgroups_factor <- factor(x = col_cellgroups_vec, levels = c("Nephron_Epithelium", "Stroma", "Immune"))


# plot heatmap body -----------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = col_fun, na_col = color_na, border = "black",
                             ## row
                             show_row_names = T, row_names_side = "left",
                             show_row_dend = F, cluster_rows = F,
                             ## column
                             show_column_names = T, column_names_side = "top",
                             column_names_gp = gpar(fontsize = 15, fontface = "bold"),
                             show_column_dend = F, 
                             column_split = col_cellgroups_factor, column_title = NULL,
                             show_heatmap_legend = F)
p
## save heatmap as png
png(filename = paste0(dir_out, "heatmap", ".png"), 
    width = 600, height = 2000, res = 150)
draw(object = p)
dev.off()


