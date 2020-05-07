# Yige Wu @WashU Apr 2020

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

# input denpendencies -----------------------------------------------------
## input log fold change
min.pct.wilcox <- 0.1
logfc.threshold.wilcox <- 0.25
logfc_by_manualsubcluster_df <- fread(input = paste0("./Resources/Analysis_Results/findmarkers/tumor_subclusters/findallmarkers_wilcox_tumorcells_by_manualsubcluster/20200427.v1/Tumormanualsubcluster.FindAllMarkers.Wilcox.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox, ".tsv"), data.table = F)
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200413.v1/meta_data.20200413.v1.tsv", data.table = F)

# set plot gene list ------------------------------------------------------
# genes_plot <- c(pbaf_genes, "SETD2", "BAP1", "KDM5C")

# make matrix for heatmap body --------------------------------------------
## add manual subcluster id
logfc_by_manualsubcluster_df$id_aliquot_wu <- mapvalues(x = logfc_by_manualsubcluster_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
logfc_by_manualsubcluster_df <- logfc_by_manualsubcluster_df %>%
  filter(p_val_adj < 0.05) %>%
  mutate(id_manualsubcluster = paste0(id_aliquot_wu, "_C", cluster))
unique(logfc_by_manualsubcluster_df$id_manualsubcluster)
plot_data_df <- dcast(data = logfc_by_manualsubcluster_df, formula = id_manualsubcluster ~ gene, value.var = "avg_logFC")
plot_data_df[is.na(plot_data_df)] <- 0
## filter column names
data_col_names <- colnames(plot_data_df)[-1]
data_col_names
data_col_names.keep <- data_col_names
## get the gene names for HIF targets
data_row_names <- plot_data_df[,1]
# data_row_names.keep <- intersect(data_row_names, genes_plot)
data_row_names.keep <- data_row_names
data_row_names.keep
## filter cell group down to epithelial cells
## reformat data frame to matrix
plot_data_mat <- as.matrix(plot_data_df[,data_col_names.keep])
plot_data_mat %>% head()
## add row names
rownames(plot_data_mat) <- data_row_names
plot_data_mat %>% head()
## filter row names
plot_data_mat <- plot_data_mat[data_row_names.keep,]
### get aliquot ids and case ids
ids_manualsubcluster <- data_row_names.keep
# ids_aliquot_wu <- str_split_fixed(string = ids_manualsubcluster, pattern = "_", n = 2)[,1]
# case_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

# make heatmap body ------------------------------------------------------------
min(plot_data_mat, na.rm = T)
max(plot_data_mat, na.rm = T)

heatmapbody_color_fun <- colorRamp2(c(-2,
                                      0,
                                      2),
                                    c("blue","white", "red"))

p <- Heatmap(matrix = plot_data_mat,
             col = heatmapbody_color_fun, column_names_gp = gpar(fontsize = 4),
             # right_annotation = row_anno, 
             # bottom_annotation = bottom_col_anno,
             # top_annotation = top_col_anno,
             # column_labels = ids_aliquot_wu,
             show_heatmap_legend = T)
p

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "tumormanualsubcluster_degs.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox, ".run", run_id, ".png")
png(filename = file2write, width = 3000, height = 1200, res = 150)
print(p)
dev.off()

