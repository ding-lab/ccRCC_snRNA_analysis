# Yige Wu @WashU May 2022

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
  "ggrastr",
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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input -------------------------------------------------------------------
cp_shRNA_deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/DEG_analysis/run_edgeR_DEG_2caki1_cp_vs_2caki1_NT/20220517.v1/Caki1_CP_vs_Caki1_NT.DEGs.20220517.v1.tsv")
# tumor_vs_pt_degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_LR_all_ccRCC_vs_pt_on_katmai/20210824.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")
tumor_vs_pt_degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/overlap_tumor_vs_normal_snRNA_bulkRNA_protein_DEGs/20210511.v1/Tumor_vs_PT_DEGs.Overlap.snRNA.bulkRNA.Protein.20210511.v1.tsv")

# filter for genes down-regulated by CP shRNA, but up in tumor cells vs. pt cells ------------------------------------------------------------------
genes_tumorcells_up <- tumor_vs_pt_degs_df$genesymbol_deg[tumor_vs_pt_degs_df$p_val_adj < 0.05 & tumor_vs_pt_degs_df$avg_log2FC > 0]
genes_tumorcells_up <- tumor_vs_pt_degs_df$genesymbol_deg[tumor_vs_pt_degs_df$direction.snRNA == "Up"]
genes_tumorcells_up
cp_shRNA_down_deg_filtered_df <- cp_shRNA_deg_df %>%
  filter(logFC < 0 ) %>%
  filter(external_gene_name %in% genes_tumorcells_up)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CP_shRNA_down_degs.potential_pathogenic.", run_id, ".tsv")
write.table(x = cp_shRNA_down_deg_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
