# Yige Wu @WashU Jun 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ggpubr",
  "ggrepel"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input -------------------------------------------------------------------
tumorcells_vs_n_nontumor_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/findmarker_wilcox_tumortissue_tumorcells_vs_normaltissue_normalcells_34samples_katmai/20220624.v1/logfcthreshold.0.minpct.0.mindiffpct.0.tsv")
tumormarkers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/tumor_specific_markers/20210701.v1/ccRCC_cells_specific_DEG_with_surface_annotations_from_3DB.txt")

# filter ------------------------------------------------------------------
genes_filter <- tumorcells_vs_n_nontumor_df$gene_symbol[tumorcells_vs_n_nontumor_df$p_val_adj < 0.05 & tumorcells_vs_n_nontumor_df$avg_log2FC > 0]
tumormarkers_filtered_df <- tumormarkers_df %>%
  filter(Gene %in% genes_filter)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_cells_specific_DEG_with_surface_annotations_from_3DB.filtred.txt")
write.table(x = tumormarkers_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
