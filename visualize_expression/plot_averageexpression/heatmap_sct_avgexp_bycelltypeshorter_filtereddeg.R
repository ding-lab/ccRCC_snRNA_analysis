# Yige Wu @WashU Aug 2020

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
## input the variable gene list
genes_df <- fread(input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/filter_markers/filter_markers_wilcox_bygroup/20200904.v1/findallmarkers_wilcox_bycellgroup.pos.logfcthreshold0.1.minpct0.1.mindiffpct0.1.top50avg_logFC.tsv", data.table = F)
## input the average expression calculated (RNA)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_bycelltypeshorter_on_katmai/20200904.v1/averageexpression_SCT_bycelltype.shorter.31_aliquot_integration.20200904.v1.tsv", data.table = F)
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200904.v1/31Aliquot.Barcode2CellType.20200904.v1.tsv", data.table = F)
## input gene annotation
genes_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_ccrcc_pathogenic_pathway_genes/20200907.v1/ccRCC_Pathogenic_Pathways_Genes.20200907.v1.tsv")

# specify the genes to show -----------------------------------------------
genes2filter <- genes_df$row_name
genes2filter <- c(genes2filter, ccRCC_drivers)

# calculate the cell counts and specify the cell types to show ------------
cellcount_df <- barcode2celltype_df %>%
  select(Cell_type.shorter) %>%
  table() %>%
  as.data.frame() %>%
  rename(Cell_type.shorter = '.') %>%
  mutate(Keep = (Freq >= 50 & Cell_type.shorter != "Unknown")) %>%
  mutate(colname_celltype = gsub(x = Cell_type.shorter, pattern = " |\\+|\\/|\\-", replacement = "."))
data_col_names.keep <- cellcount_df$colname_celltype[cellcount_df$Keep]
data_col_names.keep
# format the column names to only aliquot id ------------------------------
## filtr the rows
plot_data_df <- avgexp_df %>%
  rename(gene = V1) %>%
  filter(gene %in% genes2filter)
## remove teh prefix from the column names
data_col_names <- colnames(plot_data_df)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
data_col_names.changed
## rename the data frame
colnames(plot_data_df) <- c("gene", data_col_names.changed)
## filter the columns and make data matrix
plot_data_mat <- t(as.matrix(plot_data_df[, data_col_names.keep]))
## add gene name
colnames(plot_data_mat) <- plot_data_df$gene
plot_data_mat %>% head()

# get dimension names -----------------------------------------------------
genes_plot <- colnames(plot_data_mat)
genes_plot
# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color for the pathway annotation
colors_vhlhif <- c("TRUE" = "black", "FALSE" = "white")
## make color function for heatmap body colors
col_fun = colorRamp2(c(-1.5, 
                       0, 
                       1.5), 
                     c("blue", "white", "red"))
## make colors for the discrete ranges
colors_numbercellrange_vec <- RColorBrewer::brewer.pal(n = 6, name = "PuBuGn")
names(colors_numbercellrange_vec) <- sapply(X = seq(from = 1, to = 101, by = 20), FUN = number2rangetext)

# make annotation for cell types ------------------------------------------


# make annotation for the genes -------------------------------------------
vec_is_vhlhif <- as.character(genes_plot %in% genes_anno_df$target_genesymbol[genes_anno_df$pathway_name == "VHL-HIF"])
vec_is_epigentic_smg_related <- as.character(genes_plot %in% genes_anno_df$target_genesymbol[genes_anno_df$pathway_name == "Epigenetic machinary"])
colanno_obj = HeatmapAnnotation(VHL_HIF_Pathway = anno_simple(x = vec_is_vhlhif, col = colors_vhlhif),
                                Epigenetic_machinary_related = anno_simple(x = vec_is_epigentic_smg_related, col = colors_vhlhif))

# Heatmap -----------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
             # column_labels = ids_aliquot_wu,
             # row_labels = ids_aliquot_wu,
             # right_annotation = row_anno,
             # show_row_names = F, show_column_names = F, 
             # row_split = row_split_factor, cluster_row_slices = F, row_order = row_order_vec,
             # show_row_dend = F, row_title_rot = 0, row_title_side = "right", row_title_gp = gpar(fontsize = 25, fontface = "bold"),
             # row_gap = unit(0, "mm"),
             bottom_annotation = colanno_obj,
             # column_split = column_split_factor, cluster_column_slices = F, column_order = column_order_vec,
             # show_column_dend = F, column_title_side = "top", column_title_rot = 90, column_title_gp = gpar(fontsize = 25, fontface = "bold"),
             # column_gap = unit(0, "mm"),
             # border = "grey50",
             col = col_fun, na_col = color_na,
             show_row_names = T, 
             show_column_names = T, column_names_gp = gpar(fontsize = 5),
             show_heatmap_legend = T)
p
## save heatmap as png
png(filename = paste0(dir_out, "heatmap", ".png"), 
    width = 2000, height = 1000, res = 150)
draw(object = p)
dev.off()


