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
version_tmp <- 2
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
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200427.v1/meta_data.20200427.v1.tsv", data.table = F)
## input phenotype markers
markers_ith_df <- fread("./Resources/Analysis_Results/dependencies/write_markers_by_intratumorheterogeneity_types/20200504.v1/markergenes_by_intratumorheterogeneity_types.20200504.v1.tsv", data.table = F)

# set plot gene list ------------------------------------------------------
logfc_by_manualsubcluster_df$id_aliquot_wu <- mapvalues(x = logfc_by_manualsubcluster_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## add manual subcluster id
logfc_by_manualsubcluster_df <- logfc_by_manualsubcluster_df %>%
  mutate(id_manualsubcluster = paste0(id_aliquot_wu, "_C", cluster))
logfc_by_manualsubcluster_sig_df <- logfc_by_manualsubcluster_df %>%
  filter(gene %in% markers_ith_df$gene_symbol) %>%
  filter(p_val_adj < 0.05) %>%
  filter(avg_logFC > 0)
## get subclusters to show
ids_manualsubcluster_sig <- unique(logfc_by_manualsubcluster_sig_df$id_manualsubcluster)
ids_aliquot_wu_sig <- unique(logfc_by_manualsubcluster_df$id_aliquot_wu[logfc_by_manualsubcluster_df$id_manualsubcluster %in% ids_manualsubcluster_sig])
ids_manualsubcluster_keep <- unique(logfc_by_manualsubcluster_df$id_manualsubcluster[logfc_by_manualsubcluster_df$id_aliquot_wu %in% ids_aliquot_wu_sig])
  
# genes_plot <- c("CD44", "CXCR4", "PROM1")
genes_plot <- unique(logfc_by_manualsubcluster_sig_df$gene)

# make matrix for heatmap body --------------------------------------------
## transform
plot_data_df <- dcast(data = logfc_by_manualsubcluster_df, formula = id_manualsubcluster ~ gene, value.var = "avg_logFC")
plot_data_df[is.na(plot_data_df)] <- 0
## filter column names
data_col_names <- colnames(plot_data_df)[-1]
data_col_names
data_row_names <- plot_data_df[,1]
## filter columns
# data_col_names.keep <- data_col_names
data_col_names.keep <- genes_plot
## reformat data frame to matrix
plot_data_mat <- as.matrix(plot_data_df[,data_col_names.keep])
plot_data_mat %>% head()
## add row names
rownames(plot_data_mat) <- data_row_names
plot_data_mat %>% head()
## filter row names
data_row_names.keep <- ids_manualsubcluster_keep
data_row_names.keep
plot_data_mat <- plot_data_mat[data_row_names.keep,]
### get aliquot ids and case ids
ids_manualsubcluster <- data_row_names.keep
ids_aliquot_wu <- str_split_fixed(string = ids_manualsubcluster, pattern = "_", n = 2)[,1]
# case_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

# make column split -------------------------------------------------------
data_col_names.keep
columnsplit_vec <- mapvalues(x = data_col_names.keep, from = markers_ith_df$gene_symbol, to = as.vector(markers_ith_df$phenotype))
columnsplit_vec

# make column annotation --------------------------------------------------
genefunction_vec <- mapvalues(x = data_col_names.keep, from = markers_ith_df$gene_symbol, to = as.vector(markers_ith_df$gene_function))
regulateprocess_vec <- mapvalues(x = data_col_names.keep, from = markers_ith_df$gene_symbol, to = as.vector(markers_ith_df$regulate_process))
colors_relatedprocess_uniq <- Polychrome::dark.colors(n = length(unique(regulateprocess_vec)))
colors_relatedprocess <- mapvalues(x = regulateprocess_vec, from = unique(regulateprocess_vec), to = colors_relatedprocess_uniq)
colors_relatedprocess[is.na(regulateprocess_vec)] <- NA
col_anno <- HeatmapAnnotation(Regulate_Process = anno_text(x = regulateprocess_vec, 
                                                           location = 0.5, just = "center",
                                                    gp = gpar(fill = colors_relatedprocess, 
                                                              col = "white", 
                                                              border = "black", 
                                                              fontsize = 12, face = "bold")))

# make heatmap body ------------------------------------------------------------
min(plot_data_mat, na.rm = T)
max(plot_data_mat, na.rm = T)

heatmapbody_color_fun <- colorRamp2(c(-log(2), 0, log(2)),
                                    c("blue","white", "red"))

p <- Heatmap(matrix = plot_data_mat,
             col = heatmapbody_color_fun,
             na_col = "grey50",
             column_names_gp = gpar(fontsize = 8), 
             cluster_columns = T, column_split = columnsplit_vec, 
             # column_title_gp = gpar(fill = "white"),
             cluster_rows = F, row_split = ids_aliquot_wu,  row_title_rot = 0, 
             # row_title_gp = gpar(fill = "white"),
             border = TRUE,
             # right_annotation = row_anno, 
             bottom_annotation = col_anno,
             # top_annotation = top_col_anno,
             # column_labels = ids_aliquot_wu,
             show_heatmap_legend = T)
p

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "tumormanualsubcluster_degs.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox, ".run", run_id, ".png")
png(filename = file2write, width = 2000, height = 1200, res = 150)
print(p)
dev.off()

