# Yige Wu @WashU June 2022

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

# process ----------------------------------------------------------------
inflamm_markers_df <- inflamm_degs_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  filter(genesymbol_deg %in% tumorcell_markers_df$Gene) 
noninflamm_markers_df <- inflamm_degs_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
  filter(genesymbol_deg %in% tumorcell_markers_df$Gene) 

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "inflammatory_tumorcell_markers.", run_id, ".tsv")
write.table(x = inflamm_markers_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "noninflammatory_tumorcell_markers.", run_id, ".tsv")
write.table(x = noninflamm_markers_df, file = file2write, sep = "\t", row.names = F, quote = F)
