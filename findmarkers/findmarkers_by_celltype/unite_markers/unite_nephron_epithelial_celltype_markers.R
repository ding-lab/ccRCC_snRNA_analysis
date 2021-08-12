# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input markers for other nephron epithelial markers
dir_degs <- "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/run_bycelltype_bysample/"
deg_united_df <- NULL
for (celltype in c("Proximal tubule", "Distal convoluted tubule", "Intercalated cells", "Loop of Henle", "Principle cells", "Podocytes")) {
  celltype_print <- gsub(x = celltype, pattern = " ", replacement = "_")
  path_deg <- paste0(dir_degs, celltype_print, "/", celltype_print, "_specific_DEG_with_surface_annotations_from_3DB.txt")
  deg_df <- fread(data.table = F, input = path_deg)
  deg_df$cell_type <- celltype
  deg_united_df <- rbind(deg_united_df,deg_df)
}
table(deg_united_df$cell_type)
## input tumor-cell markers
deg_tumorcell_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/tumor_specific_markers/20210701.v1/Tumor cells_specific_DEG_with_surface_annotations_from_3DB.txt")
deg_tumorcell_df$cell_type <- "Tumor cells"
deg_united_df <- rbind(deg_united_df, deg_tumorcell_df)
table(deg_united_df$cell_type)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC.Nephron_Epithelial_CellType_AlgorithmeticMarkers.", run_id, ".tsv")
write.table(file = file2write, x = deg_united_df, row.names = F, quote = F, sep = "\t")
