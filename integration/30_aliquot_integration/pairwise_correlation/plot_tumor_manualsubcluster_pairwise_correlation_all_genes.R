# Yige Wu @WashU March 2020
## running on local
## for calculating the aliquot-pairwise correlation coefficients for averaged expression of all genes

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/plotting.R")
library(ComplexHeatmap)
library(circlize)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input id meta data
id_metadata_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## input the spearman pairwise correlation result
pearson_coef_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/pairwise_correlation/calculate_tumor_manualsubcluster_pairwise_correlation_all_genes/20200407.v1/avg_exp_by_tumorsubluster.all_genes.pearson_coef20200407.v1.tsv", data.table = F)
## input pathway score
pathway_score_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/calculate_pathway_score/unite_pathway_scores/20200407.v1/avg_pathway_scaled_avgexp.20200407.v1.tsv", data.table = F)

# make data matrix for heatmap body ---------------------------------------
## reformat data frame to matrix
plot_data_df <- pearson_coef_df
plot_data_mat <- as.matrix(plot_data_df[,-1])
plot_data_mat %>% head()
## add row names
rownames(plot_data_mat) <- plot_data_df$V1
plot_data_mat %>% head()
### get aliquot ids and case ids
tumorsubcluster_ids <- rownames(plot_data_mat)
aliquot_ids <- str_split_fixed(string = tumorsubcluster_ids, pattern = "_", n = 2)[,1]
case_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

# make top column annotation --------------------------------------------------
## make annotation data frame from the pathway score first
top_col_anno_df <- dcast(data = pathway_score_df, formula = exp_manual_subcluster_name ~ pathway_name, value.var = "pathway_score")
rownames(top_col_anno_df) <- top_col_anno_df$exp_manual_subcluster_name
## filter by tumor subcluster id 
top_col_anno_df <- top_col_anno_df[tumorsubcluster_ids,]
### make color for hif pathway score
hifscore_color_fun <-  colorRamp2(c(quantile(top_col_anno_df$HIF, 0.1, na.rm=T), 
                                    quantile(top_col_anno_df$HIF, 0.5, na.rm=T), 
                                    quantile(top_col_anno_df$HIF, 0.9, na.rm=T)),
                                  c("blue", "white", "red"))
### make color for PBAF complex score
pbafscore_color_fun <-  colorRamp2(c(quantile(top_col_anno_df$PBAF, 0.1, na.rm=T), 
                                     quantile(top_col_anno_df$PBAF, 0.5, na.rm=T), 
                                     quantile(top_col_anno_df$PBAF, 0.9, na.rm=T)),
                                   c("blue", "white", "red"))
### make color for mTOR patway score
mtorscore_color_fun <-  colorRamp2(c(quantile(top_col_anno_df$MTOR, 0.1, na.rm=T), 
                                     quantile(top_col_anno_df$MTOR, 0.5, na.rm=T), 
                                     quantile(top_col_anno_df$MTOR, 0.9, na.rm=T)),
                                   c("blue", "white", "red"))
### make top column annotation object
top_col_anno = HeatmapAnnotation(HIF_Pathway_Score = anno_simple(x = top_col_anno_df$HIF,
                                                                 simple_anno_size = unit(3, "mm"),
                                                                 col = hifscore_color_fun),
                                 PBAF_Complex_Score = anno_simple(x = top_col_anno_df$PBAF,
                                                                  simple_anno_size = unit(3, "mm"),
                                                                  col = pbafscore_color_fun),
                                 mTOR_Pathway_Score = anno_simple(x = top_col_anno_df$MTOR,
                                                                  simple_anno_size = unit(3, "mm"),
                                                                  col = mtorscore_color_fun))

# make row annotation -------------------------------------------
### create row annotation
row_anno = rowAnnotation(foo = anno_text(case_ids, 
                                         location = 0.5, just = "center",
                                         gp = gpar(fill = uniq_case_colors[case_ids], col = "white", border = "black"),
                                         width = max_text_width(case_ids)*1.2))


# make bottom column annotation -------------------------------------------
bottom_col_anno = HeatmapAnnotation(foo = anno_text(case_ids, 
                                                    location = 0.5, just = "center",
                                                    gp = gpar(fill = uniq_case_colors[case_ids], col = "white", border = "black"),
                                                    width = max_text_width(case_ids)*1.2))



# plot pearson pairwise correlation for variably expressed genes within tumor cells ------------------------------------------------------
## make color function for heatmap body colors
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
## make heatmap
p <- Heatmap(matrix = plot_data_mat,
             column_km = 2, column_km_repeats = 100,
             row_km = 2, row_km_repeats = 100,
             col = col_fun, 
             right_annotation = row_anno,
             bottom_annotation = bottom_col_anno, 
             top_annotation = top_col_anno,
             show_heatmap_legend = F)
p
## make legend for heattmap body
heatmap_lgd = Legend(col_fun = col_fun, title = "Pearson's coeffcient\n(variably expressed genes\nwithin tumor cells)", direction = "vertical")
## make legend for top annotation
annotation_lgd = list(
  heatmap_lgd,
  Legend(labels = c("copy loss", "copy gain"), 
         title = "Chr Arm Copy Number Alteration", 
         legend_gp = gpar(fill = c("blue", "red"), direction = "horizontal"),
         direction = "horizontal"),
  Legend(labels = names(variant_class_colors), 
         title = "Mutation Class", 
         legend_gp = gpar(fill = variant_class_colors),
         direction = "horizontal"))

## save heatmap
png(filename = paste0(dir_out, "avg_exp.all_genes.pearson_coef.heatmap.", run_id, ".png"), 
    width = 4500, height = 4000, res = 150)
### combine heatmap and heatmap legend
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()
