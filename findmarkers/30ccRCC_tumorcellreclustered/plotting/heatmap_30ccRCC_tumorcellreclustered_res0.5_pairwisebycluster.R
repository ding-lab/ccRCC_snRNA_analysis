# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "circlize",
  "ComplexHeatmap"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input DEGs
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/30ccRCC_tumorcellreclustered/findmarker_30ccRCC_tumorcellreclustered_res0/20220405.v1/tumorcellsreclustered.pairwisebycluster.markers.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")
## input the barcode to cluster id
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_30_ccRCC_tumorcells/20220404.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.20220404.v1.tsv")

# make matrix data -----------------------------------------------------------------
plotdata_df <- results_df %>%
  filter(p_val_adj < 0.05) %>%
  mutate(diff_pct = (pct.1 - pct.2)) %>%
  filter(diff_pct >= 0.1) %>%
  mutate(x_plot = ifelse(ident.1 > ident.2, ident.1, ident.2)) %>%
  mutate(y_plot = ifelse(ident.1 > ident.2, ident.2, ident.1)) %>%
  group_by(x_plot, y_plot) %>%
  summarise(number_degs = n())

plotdata_df <- rbind(plotdata_df,
                     data.frame(x_plot = c(1, 0), y_plot = c(1, 0), number_degs = c(0, 0)))
plotdata_df$x_plot <- factor(plotdata_df$x_plot)
plotdata_df$y_plot <- factor(plotdata_df$y_plot)
plotdata_wide_df <- dcast(data = plotdata_df, formula = x_plot ~ y_plot, value.var = "number_degs")
plotdata_mat <- as.matrix(plotdata_wide_df[,-1])
rownames(plotdata_mat) <- plotdata_wide_df$x_plot

# make colors -------------------------------------------------------------
col_fun = circlize::colorRamp2(breaks = c(0, 50, 500), colors = c("blue", "white", "red"))

# make column annotation -----------------------------------------------------
anno_df <- barcode2cluster_df %>%
  group_by(ident) %>%
  summarise(number_cells = n())
anno_obj <- rowAnnotation(Numberof_cells = anno_df$number_cells)

# plot --------------------------------------------------------------------
Heatmap(matrix = plotdata_mat, col = col_fun, cluster_columns = F, cluster_rows = F, left_annotation = anno_obj, row_names_side = "left")
