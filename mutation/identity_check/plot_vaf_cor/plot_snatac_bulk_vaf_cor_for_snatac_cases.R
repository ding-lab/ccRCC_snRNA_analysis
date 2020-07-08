# Yige Wu @WashU Jul 2020

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

# input dependencies -----------------------------------------------------
## input correlation coefficient
plot_data_mat <- readRDS(file = "./Resources/Analysis_Results/mutation/identity_check/calculate_vaf_cor/calculate_snatac_bulk_vaf_cor_for_snatac_samples/20200706.v1/Correlation_Coeffcients.snATAC_vs_bulk.snATAC_Cases.20200706.v1.RDS")

# plot heatmap body with white-yellow-red ------------------------------------------------------
## make color function for heatmap body colors
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
## make heatmap
p <- Heatmap(matrix = plot_data_mat, 
             # show_row_names = F, show_column_names = F, show_column_dend = F, show_row_dend = F,
             # column_km = 2, column_km_repeats = 100,
             # row_km = 2, row_km_repeats = 100,
             # right_annotation = row_anno, 
             # bottom_annotation = bottom_col_anno, 
             # top_annotation = top_col_anno,
             col = col_fun, name = "Correlation Coeffcient",
             show_heatmap_legend = T)
p

# write output ------------------------------------------------------------
## save heatmap as png
file2write <- paste0(dir_out, "snatac_bulk_vaf_cor_for_snatac_samples.", ".png")
png(filename = file2write, 
    width = 1000, height = 800, res = 150)
print(p)
dev.off()


