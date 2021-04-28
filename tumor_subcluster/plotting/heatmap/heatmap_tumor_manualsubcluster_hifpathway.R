# Yige Wu @WashU Apr 2020
## running on local
## for plotting average expression of known pathogenic pathway genes for each tumor subclusters (manually grouped)
## VHL-HIF pathway

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input scaled average expression
scaled_avgexp_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/other/scale_averageexpression/20200407.v1/averageexpression_tumor_cells_by_manual_subcluster.scaled.20200407.v1.tsv", data.table = F)
## input pathway members
ora_pathwaymember_list <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/calculate_pathway_score/calculate_ora_pathway_score/20200407.v1/ora_pathway.pathway_members.20200407.v1.RDS")
## input id meta data
id_metadata_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## input cnv fraction by chr
frac_cnv_wide_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/summarize_cnv_fraction/estimate_fraction_of_tumorcells_with_expectedcnv_perchrregion_per_manualsubcluster_using_cnvgenes/20200407.v1/fraction_of_tumorcells.expectedCNA.by_chr_region.by_manual_subcluster.20200407.v1.tsv", data.table = F)
rownames(frac_cnv_wide_df) <- frac_cnv_wide_df$manual_cluster_name

# make data matrix for heatmap body ---------------------------------------
## format the column names to only aliquot id 
data_col_names <- colnames(scaled_avgexp_df)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
## rename the data frame
colnames(scaled_avgexp_df) <- c("gene", data_col_names.changed)
## filter out unwanted manual clusters
data_col_names.keep <- data_col_names.changed[!grepl(pattern = "MCNA", x = data_col_names.changed)]
## reformat data frame to matrix
plot_data_df <- scaled_avgexp_df
plot_data_mat <- as.matrix(plot_data_df[,data_col_names.keep])
plot_data_mat %>% head()
## add row names
rownames(plot_data_mat) <- plot_data_df$gene
plot_data_mat %>% head()
### get aliquot ids and case ids
tumorsubcluster_ids <- colnames(plot_data_mat)
aliquot_ids <- str_split_fixed(string = tumorsubcluster_ids, pattern = "_", n = 2)[,1]
case_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

# make top column annotation --------------------------------------------------
## make annotation data frame with copy number profile first
top_col_anno_df_wide <- frac_cnv_wide_df
top_col_anno_df <- top_col_anno_df_wide[,-1]
rownames(top_col_anno_df) <- top_col_anno_df_wide$manual_cluster_name
top_col_anno_df <- top_col_anno_df[colnames(plot_data_mat),]
top_col_anno_df[is.na(top_col_anno_df)] <- 0

### make top column annotation object
top_col_anno = HeatmapAnnotation(Fraction_Cells_With_3p_Loss = anno_barplot(x = top_col_anno_df$`3p`, 
                                                                           gp = gpar(fill = 1, col = 1)),
                                 Fraction_Cells_With_MDM4_Gain = anno_barplot(x = top_col_anno_df$`1q32`, 
                                                                           gp = gpar(fill = 1, col = 1)),
                                 show_legend = T)

# make bottom column annotation -------------------------------------------
bottom_col_anno = HeatmapAnnotation(foo = anno_text(case_ids, 
                                                    location = 0.5, just = "center",
                                                    gp = gpar(fill = uniq_case_colors[case_ids], col = "white", border = "black"),
                                                    width = max_text_width(case_ids)*1.2))

# plot pathway genes  ------------------------------------------------------
## get the subset of genes
plot_genes_tmp <- ora_pathwaymember_list[["HIF-1-alpha transcription factor network"]]
# plot_genes_tmp <- c(plot_genes_tmp, as.vector(ccrcc_cna_genes_df$gene_symbol[ccrcc_cna_genes_df$chr_region == "3p"]))
# plot_genes_tmp <- intersect(plot_genes_tmp, variable_genes_df$gene)
## get the subset of plot data
plot_data_mat_tmp <- plot_data_mat[plot_genes_tmp,]
## make function for colors
heatmapbody_color_fun <- colorRamp2(c(quantile(plot_data_mat_tmp, 0.1, na.rm=T), 
                                      quantile(plot_data_mat_tmp, 0.5, na.rm=T), 
                                      quantile(plot_data_mat_tmp, 0.9, na.rm=T)),
                                    c("blue", "white", "red"))

p <- Heatmap(matrix = plot_data_mat_tmp,
             col = heatmapbody_color_fun,
             bottom_annotation = bottom_col_anno,
             top_annotation = top_col_anno,
             show_heatmap_legend = T)
p
## save heatmap
png(filename = paste0(dir_out, "avg_exp.hif_pathway.by_tumor_manualsubcluster.heatmap.", run_id, ".png"), 
    width = 4500, height = 1500, res = 150)
### combine heatmap and heatmap legend
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()


