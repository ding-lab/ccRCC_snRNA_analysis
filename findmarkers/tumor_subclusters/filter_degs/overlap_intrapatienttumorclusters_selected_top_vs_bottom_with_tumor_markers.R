# Yige Wu @WashU Dec 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "data.table",
  "stringr"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input data --------------------------------------------------------------
inflamm_degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_intrapatienttumorclusters_inflammatory_top_vs_bottom_katmai/20220610.v1/Inflammatory_score_top_vs_bottom_tumorclusters..logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")
## input tumor0cell markers
tumorcell_markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/tumor_specific_markers/20210701.v1/ccRCC_cells_specific_DEG_with_surface_annotations_from_3DB.txt")

# input all degs ----------------------------------------------------------
dir_input <- "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_intrapatienttumorclusters_selectedpathway_top_vs_bottom_katmai/20221206.v1/"
files_input <- list.files(path = dir_input)
degs_all_df <- NULL
for (file_tmp in files_input) {
  degs_tmp <- fread(data.table = F, input = paste0(dir_input, file_tmp))
  geneset_tmp <- str_split(string = file_tmp, pattern = "\\.")[[1]][1]
  degs_tmp$geneset <- geneset_tmp
  degs_all_df <- rbind(degs_all_df, degs_tmp)
}

# overlap ----------------------------------------------------------------
degs_overlap_df <- degs_all_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  filter(genesymbol_deg %in% tumorcell_markers_df$Gene) 
table(degs_overlap_df$geneset)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "selectedpathway_markers_overlap_tumorcell_markers.", run_id, ".tsv")
write.table(x = degs_overlap_df, file = file2write, sep = "\t", row.names = F, quote = F)
