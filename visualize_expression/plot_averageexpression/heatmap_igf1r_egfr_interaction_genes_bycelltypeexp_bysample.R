# Yige Wu @WashU OCt 2020

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
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_byallacelltype_byaliquot_on_katmai/20201022.v1/avgexp.SCT.scale.data.bycelltypeshorter.byaliquot.31_aliquot_integration.20201022.v1.tsv", data.table = F)
## input interaction strength
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/other/scale_pair_cell_types_vs_sample_sig_mean/20200924.v1/cellphonedb.pair_cell.types_vs_sample.sig_mean.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# specify the genes to show -----------------------------------------------
celltype1 <- "Macrophages"
celltype2 <- "Tumor.cells"
genes_celltype1 <- c("IGF1", "HBEGF")
genes_celltype2 <-c("IGF1R", "EGFR")
genes2filter <- c(genes_celltype1, genes_celltype2)
rownames_sorted <- c("IGF1_Macrophages", "IGF1R_Tumor.cells",
                     "HBEGF_Macrophages", "EGFR_Tumor.cells")
# get interaction strength and sorted easy ids -----------------------------------------------------
strength_filtered_df <- cellphone_df %>%
  filter(pair_cell.types %in% c("IGF1_IGF1R.Macrophages|Tumor cells", "EGFR_HBEGF.Tumor cells|Macrophages"))
strength_sort_df <- data.frame(t(strength_filtered_df[,-1]))
colnames(strength_sort_df) <- strength_filtered_df$pair_cell.types
strength_sort_df$number_sig <- rowSums(!is.na(strength_sort_df))
strength_sort_df$Easy_id <- colnames(strength_filtered_df)[-1]
strength_sort_df <- strength_sort_df %>%
  arrange(number_sig, desc(`IGF1_IGF1R.Macrophages|Tumor cells`), `EGFR_HBEGF.Tumor cells|Macrophages`)
## get sorted easy ids
easyids_sorted <- strength_sort_df$Easy_id

# format the column names to only aliquot id ------------------------------
## filtr the rows
plot_data_df <- avgexp_df %>%
  filter(V1 %in% genes2filter)
rownames(plot_data_df) <- plot_data_df$V1
## filter the columns and make data matrix
### get the column names to keep
sct_aliquot_celltype_vec <- colnames(plot_data_df)[-1]
colname_mapping_df <- data.frame(sct.aliquot_celltype = sct_aliquot_celltype_vec)
colname_mapping_df <- colname_mapping_df %>%
  mutate(aliquot_celltype = gsub(x = sct.aliquot_celltype, pattern = "SCT.", replacement = "")) %>%
  mutate(aliquot = str_split_fixed(string = aliquot_celltype, pattern = "_", n = 2)[,1]) %>%
  mutate(celltype = str_split_fixed(string = aliquot_celltype, pattern = "_", n = 2)[,2]) %>%
  arrange(celltype, aliquot)
rownames(colname_mapping_df) <- colname_mapping_df$sct.aliquot_celltype
### get ordered column names by cell type  
colnames_celltype1 <- colname_mapping_df$sct.aliquot_celltype[colname_mapping_df$celltype == celltype1];colnames_celltype1 <- as.vector(colnames_celltype1); 
colnames_celltype2 <- colname_mapping_df$sct.aliquot_celltype[colname_mapping_df$celltype == celltype2];colnames_celltype2 <- as.vector(colnames_celltype2); colnames_celltype2
### get different ids
aliquot_1diff2 <- setdiff(colname_mapping_df$aliquot[colname_mapping_df$celltype == celltype1], 
                          colname_mapping_df$aliquot[colname_mapping_df$celltype == celltype2])
aliquot_2diff1 <- setdiff(colname_mapping_df$aliquot[colname_mapping_df$celltype == celltype2], 
                          colname_mapping_df$aliquot[colname_mapping_df$celltype == celltype1])
plot_data_df1 <- plot_data_df[genes_celltype1, colnames_celltype1]; colnames(plot_data_df1) <- colname_mapping_df[colnames_celltype1, "aliquot"]
plot_data_df2 <- plot_data_df[genes_celltype2, colnames_celltype2]; colnames(plot_data_df2) <- colname_mapping_df[colnames_celltype2, "aliquot"]

if (length(aliquot_1diff2) > 0) {
  plot_data_df2[, aliquot_1diff2] <- NA
}
if (length(aliquot_2diff1) > 0) {
  plot_data_df1[, aliquot_2diff1] <- NA
}
plot_data_mat <- as.matrix(rbind(plot_data_df1,
                                 plot_data_df2))
## add gene name
rownames(plot_data_mat) <- c(paste0(genes_celltype1, "_", celltype1),
                             paste0(genes_celltype2, "_", celltype2))
### sort rows
plot_data_mat <- plot_data_mat[rownames_sorted,]
## get dimension names
easyids_col <- mapvalues(x = colnames(plot_data_mat), from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## rename column names
colnames(plot_data_mat) <- easyids_col
### sort columns
plot_data_mat <- plot_data_mat[, easyids_sorted[easyids_sorted %in% easyids_col]]

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plot_data_mat))
colors_heatmapbody = colorRamp2(c(-1.5, 
                       0, 
                       1.5), 
                     c(color_blue, "white", color_red))
colors_interaction <- colorRamp2(c(0, 2, 4, 6, 8, 10), 
                                 c(brewer.pal(n = 6, name = "YlOrRd")))

# make column annotatin ---------------------------------------------------
igf_igf1r_interaction_vec <- strength_sort_df$`IGF1_IGF1R.Macrophages|Tumor cells`[strength_sort_df$Easy_id %in% easyids_col]
hbegf_egfr_interaction_vec <- strength_sort_df$`EGFR_HBEGF.Tumor cells|Macrophages`[strength_sort_df$Easy_id %in% easyids_col]

colanno_obj = HeatmapAnnotation(
  IGF1_IGF1R_Interaction = anno_simple(x = igf_igf1r_interaction_vec, col = colors_interaction),
  HBEGF_EGFR_Interaction = anno_simple(x = hbegf_egfr_interaction_vec, col = colors_interaction),
  annotation_name_gp = gpar(fontsize = 15, fontface = "italic"), annotation_name_side = "left")


# Heatmap -----------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody, na_col = color_na, border = "black",
                             ## column
                             cluster_columns = F, 
                             top_annotation = colanno_obj,
                             ## row
                             cluster_rows = F, row_names_side = "left",
                             show_heatmap_legend = F)
p

# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_interaction, 
         title = "Interaction strength", 
         title_gp = gpar(fontsize = 10, fontface = "bold"),
         legend_width = unit(4, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_heatmapbody, 
         title = "Average expression", 
         title_gp = gpar(fontsize = 10, fontface = "bold"),
         legend_width = unit(4, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "interactions", ".png")
png(file2write, width = 1500, height = 500, res = 150)
draw(object = p, 
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()

file2write <- paste0(dir_out, "interactions", ".pdf")
pdf(file2write, width = 10, height = 3.75, useDingbats = F)
draw(object = p, 
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()

