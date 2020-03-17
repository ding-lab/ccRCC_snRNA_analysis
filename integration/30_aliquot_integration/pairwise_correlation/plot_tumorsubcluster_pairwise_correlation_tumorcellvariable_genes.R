# Yige Wu @WashU March 2020
## running on local
## for calculating the aliquot-pairwise correlation coefficients for averaged expression of all genes

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/aes.R")
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
## input cnv profile
bulk_cnv_profile_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/bulk/write_sample_bicseq_cnv_profile/20200227.v1/Bulk_WGS_Chr_CNV_Profile.20200227.v1.tsv", data.table = F)
## input the spearman pairwise correlation result
pearson_coef.tumorcellvariable_genes.df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/pairwise_correlation/calculate_tumorsubcluster_pairwise_correlation_tumorcellvariable_genes/20200312.v1/avg_exp_by_tumorsubluster.tumorcellvaraible_genes.pearson_coef20200312.v1.tsv", data.table = F)
## input mutation
maf_df <- loadMaf()

# make data matrix for heatmap body ---------------------------------------
## reformat data frame to matrix
plot_data_df <- pearson_coef.tumorcellvariable_genes.df
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
## make annotation data frame with copy number profile first
top_col_anno_df <- bulk_cnv_profile_df %>%
  select('3p', '5q', "14q")
rownames(top_col_anno_df) <- bulk_cnv_profile_df$Case
top_col_anno_df <- top_col_anno_df[case_ids,]
rownames(top_col_anno_df) <- rownames(plot_data_mat)
## make annotation data frame with mutation type
### filter maf file by case id
maf_df <- maf_df %>%
  mutate(case_id = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", n = 2)[,1]) %>%
  filter(case_id %in% case_ids)
var_class_df <- generate_somatic_mutation_matrix(pair_tab = ccRCC_SMGs, maf = maf_df)
var_class_anno_df <- t(var_class_df[,-1])
var_class_anno_df
var_class_anno_df <- var_class_anno_df[case_ids,]
var_class_anno_df
rownames(var_class_anno_df) <- rownames(plot_data_mat)
### combine mutation data into top column annotation
top_col_anno_df <- cbind(top_col_anno_df, var_class_anno_df)
### remove variant info for the normal sample
normal_aliquot_ids <- id_metadata_df$Aliquot.snRNA[id_metadata_df$Sample_Type == "Normal"]
top_col_anno_df[rownames(top_col_anno_df) %in% normal_aliquot_ids,] <- ""
### make color for this annotation data frame
top_col_anno_colors <- lapply(colnames(top_col_anno_df), function(g) {
  if (g %in% c('3p', '5q', "14q")) {
    color_vector <- c("gain" = "red", "loss" = "blue")
  }
  if (g %in% ccRCC_SMGs) {
    color_vector <- variant_class_colors
  }
  return(color_vector)
})
names(top_col_anno_colors) <- colnames(top_col_anno_df)
### make top column annotation object
top_col_anno = HeatmapAnnotation(df = top_col_anno_df, col = top_col_anno_colors, show_legend = F)

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
## make function for colors
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
## make heatmap
p <- Heatmap(matrix = plot_data_mat,
             col = col_fun, 
             right_annotation = row_anno, 
             bottom_annotation = bottom_col_anno, 
             # top_annotation = top_col_anno,
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
png(filename = paste0(dir_out, "avg_exp.tumorcellvariable_genes.pearson_coef.heatmap.", run_id, ".png"), 
    width = 5500, height = 5000, res = 150)
### combine heatmap and heatmap legend
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()
