# Yige Wu @WashU March 2020

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
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv", data.table = F)
## input the spearman pairwise correlation result
pearson_coef.tumorcellvariable_genes.df <- fread(input = "./Resources/Analysis_Results/integration/30_aliquot_integration/pairwise_correlation/calculate_tumor_manualsubcluster_pairwise_correlation_tumorcellvariable_genes/20200721.v1/avg_exp_by_tumorsubluster.tumorcellvaraible_genes.pearson_coef20200721.v1.tsv", data.table = F)

# make data matrix for heatmap body ---------------------------------------
## reformat data frame to matrix
plot_data_df <- pearson_coef.tumorcellvariable_genes.df
plot_data_mat <- as.matrix(plot_data_df[,-1])
plot_data_mat %>% head()
## add row names
rownames(plot_data_mat) <- plot_data_df$V1
plot_data_mat %>% head()

# get ids and filter data matrix -----------------------------------------------------------------
## get tumor subcluster names
ids_tumorsubcluster <- rownames(plot_data_mat)
names_tumorsubclusters <- gsub(x = ids_tumorsubcluster, pattern = "\\.", replacement = "-")
names_tumorsubclusters
## get case ids
ids_case <- str_split_fixed(string = names_tumorsubclusters, pattern = "-T", n = 2)[,1]
ids_case_plot_row_uniq <- c("C3L-00088")
ids_case_plot_column_uniq <- c("C3L-00088")
## filter data matrix
ids_tumorsubcluster_plot_row <- ids_tumorsubcluster[ids_case %in% ids_case_plot_row_uniq]
ids_tumorsubcluster_plot_column <- ids_tumorsubcluster[ids_case %in% ids_case_plot_column_uniq]
plot_data_mat <- plot_data_mat[ids_tumorsubcluster_plot_row, ids_tumorsubcluster_plot_column]
## get case ids to plot
ids_case_plot_row <- ids_case[ids_case %in% ids_case_plot_row_uniq]
ids_case_plot_column <- ids_case[ids_case %in% ids_case_plot_column_uniq]
## get tumor subcluster names to plot
names_tumorsubclusters_plot_row <- names_tumorsubclusters[ids_case %in% ids_case_plot_row_uniq]
names_tumorsubclusters_plot_column <- names_tumorsubclusters[ids_case %in% ids_case_plot_column_uniq]
## get readable sample id
ids_aliquot_wu_plot_row <- str_split_fixed(string = names_tumorsubclusters_plot_row, pattern = "_", n = 2)[,1]
ids_aliquot_wu_plot_column <- str_split_fixed(string = names_tumorsubclusters_plot_column, pattern = "_", n = 2)[,1]
## get suffixes of the aliquot ids
suffixes_aliquot_id_plot_row <- str_split_fixed(string = ids_aliquot_wu_plot_row, pattern = "-", n = 3)[,3]
suffixes_aliquot_id_plot_column <- str_split_fixed(string = ids_aliquot_wu_plot_column, pattern = "-", n = 3)[,3]
## get the cluster names
names_cluster_suffix_plot_row <- str_split_fixed(string = names_tumorsubclusters_plot_row, pattern = "_", n = 2)[,2]
names_cluster_suffix_plot_column <- str_split_fixed(string = names_tumorsubclusters_plot_column, pattern = "_", n = 2)[,2]

# process bulk omics ------------------------------------------------------
## filter samples with either snRNA data
id_metadata_filtered_df <- idmetadata_df %>%
  filter(snRNA_available) %>%
  filter(Case %in% unique(c(ids_case_plot_row_uniq, ids_case_plot_column_uniq)))

# sort case ids -----------------------------------------------------------
## order case ids
factor_ids_case <- factor(x = ids_case, levels = clustercount_by_case$ids_case)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make colors for histogical type
colors_hist_type <- c("Clear cell renal cell carcinoma" = "#fc8d62", "non-Clear cell renal cell carcinoma" = "#8da0cb")
## make color function for heatmap body colors
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
## make color for the cluster name suffix
unique(names_cluster_suffix)
colors_clustername <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
names(colors_clustername) <- paste0("C", 1:8)
swatch(colors_clustername)

# make row annotation -------------------------------------------
## sort omics
## create left row annotation
row_anno = rowAnnotation(Sample_Type_Suffix = anno_text(x = suffixes_aliquot_id_plot_row, 
                                                        location = 0.5, just = "center",
                                                        gp = gpar(col = "black", fontsize = 15, fontface = "bold", fill = colors_tumor_segments[suffixes_aliquot_id_plot_row]), 
                                                        width = unit(1, "cm"),
                                                        rot = 0),
                         Cluster_Name = anno_text(x = names_cluster_suffix_plot_row, 
                                                  location = 0.5, just = "center",
                                                  gp = gpar(col = "black", fontsize = 15, fontface = "bold", fill = colors_clustername[names_cluster_suffix_plot_row]), 
                                                  width = unit(1, "cm"),
                                                  rot = 0),
                         annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = 20))

# make column annotation -------------------------------------------
## do not show T1 for samples with only one segments
col_anno = HeatmapAnnotation(Sample_Type_Suffix = anno_text(x = suffixes_aliquot_id_plot_column, 
                                                            location = 0.5, just = "center",
                                                            gp = gpar(col = "black", fontsize = 15, fontface = "bold", fill = colors_tumor_segments[suffixes_aliquot_id_plot_column]), 
                                                            height = unit(1, "cm"),
                                                            rot = 0),
                             Cluster_Name = anno_text(x = names_cluster_suffix_plot_column, 
                                                      location = 0.5, just = "center",
                                                      gp = gpar(col = "black", fontsize = 15, fontface = "bold", fill = colors_clustername[names_cluster_suffix_plot_row]), 
                                                      height = unit(1, "cm"),
                                                      rot = 0))

# plot pearson pairwise correlation for variably expressed genes within tumor cells ------------------------------------------------------
## make heatmap
p <- Heatmap(matrix = plot_data_mat,
             width = unit(nrow(plot_data_mat), "cm"), height = unit(ncol(plot_data_mat), "cm"),
             col = col_fun, 
             cluster_rows = T,
             show_row_dend = F,
             row_gap = unit(0, "mm"),
             left_annotation = row_anno,
             cluster_columns = T,
             show_column_dend = F,
             column_gap = unit(0, "mm"),
             border = "grey50",
             bottom_annotation= col_anno,
             show_row_names = F,
             show_column_names = F,
             show_heatmap_legend = F)
## save heatmap
file2write <- paste0(dir_out, ids_case_plot_row_uniq, ".pdf")
pdf(file2write,
    width = 7, height = 7)
draw(object = p)
dev.off()