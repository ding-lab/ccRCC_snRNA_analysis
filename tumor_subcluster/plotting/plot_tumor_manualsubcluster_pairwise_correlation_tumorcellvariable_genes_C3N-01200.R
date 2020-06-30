# Yige Wu @WashU March 2020
## running on local
## for calculating the aliquot-pairwise correlation coefficients for averaged expression of all genes

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
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## input the spearman pairwise correlation result
pearson_coef.tumorcellvariable_genes.df <- fread(input = "./Resources/Analysis_Results/integration/30_aliquot_integration/pairwise_correlation/calculate_tumor_manualsubcluster_pairwise_correlation_tumorcellvariable_genes/20200325.v1/avg_exp_by_tumorsubluster.tumorcellvaraible_genes.pearson_coef20200325.v1.tsv", data.table = F)
## set case to do
id_case_plot <- "C3N-01200"

# rename the input matrix -------------------------------------------------
ids_tumorsubcluster_all <- colnames(plot_data_df)[-1]
ids_aliquot_all <- str_split_fixed(string = ids_tumorsubcluster_all, pattern = "_", n = 2)[,1]
ids_aliquot_all
ids_case_all <- mapvalues(x = ids_aliquot_all, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
ids_case_all
ids_aliquot_wu_all <- mapvalues(x = ids_aliquot_all, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
ids_aliquot_wu_all
ids_cluster_all <- str_split_fixed(string = ids_tumorsubcluster_all, pattern = "_MC", n = 2)[,2]
ids_cluster_all
names_clusters_all <- paste0("C", (as.numeric(ids_cluster_all)+1))
names_clusters_all

# make data matrix for heatmap body ---------------------------------------
## reformat data frame to matrix
plot_data_df <- pearson_coef.tumorcellvariable_genes.df
## get column names to filter down to
colnames_filtered_raw <- ids_tumorsubcluster_all[ids_case_all %in% id_case_plot]
colnames_filtered_raw
## filter
plot_data_mat <- plot_data_df[plot_data_df$V1 %in% colnames_filtered_raw, colnames_filtered_raw]
## rename column and row name
colnames_filtered <- paste0(ids_aliquot_wu_all[ids_case_all %in% id_case_plot], "_", names_clusters_all[ids_case_all %in% id_case_plot])
colnames_filtered
rownames(plot_data_mat) <- colnames_filtered
colnames(plot_data_mat) <- colnames_filtered
plot_data_mat <- as.matrix(plot_data_mat)
# ### get aliquot ids and case ids
# tumorsubcluster_ids <- rownames(plot_data_mat)
# aliquot_ids <- str_split_fixed(string = tumorsubcluster_ids, pattern = "_", n = 2)[,1]
# case_ids <- mapvalues(x = aliquot_ids, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))

# plot pearson pairwise correlation for variably expressed genes within tumor cells ------------------------------------------------------
## make color function for heatmap body colors
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
## make heatmap
p <- Heatmap(matrix = plot_data_mat,
             col = col_fun, 
             show_heatmap_legend = T)
p
## save heatmap
png(filename = paste0(dir_out, "avg_exp.tumorcellvariable_genes.pearson_coef.heatmap.", run_id, ".png"), 
    width = 1000, height = 1000, res = 150)
### combine heatmap and heatmap legend
draw(object = p)
dev.off()
