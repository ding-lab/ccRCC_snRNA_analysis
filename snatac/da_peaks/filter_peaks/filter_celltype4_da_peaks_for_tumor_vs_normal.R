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
## input united DA peaks
peaks2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_celltype4_da_peaks/20200915.v1/DA_peaks.chromvar.MergedObj.byCell_group4.Annotated.20200915.v1.tsv")
## input tumor vs normal DEGs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_tumor_vs_pt_on_katmai/20200903.v1/findallmarkers_wilcox_tumorcells_vs_pt.20200903.v1.tsv")

# filter DEGs -------------------------------------------------------------
degs_filtered_df <- degs_df %>%
  filter(p_val_adj < 0.05)

# filter by promoter region and DEGs --------------------------------------
peaks_filtered_df <- peaks2celltype_df %>%
  filter(avg_logFC > 0) %>%
  filter(Cell_type.filename %in% c("Tumor cells", "Normal epithelial cells")) %>%
  filter(annotation == "Promoter")
## filter DEGs
peaks_filtered_df <- peaks_filtered_df %>%
  filter((Cell_type.filename == "Tumor cells" & SYMBOL %in% degs_filtered_df$row_name[degs_filtered_df$avg_logFC > 0]) | (Cell_type.filename == "Normal epithelial cells" & SYMBOL %in% degs_filtered_df$row_name[degs_filtered_df$avg_logFC < 0]))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DA_peaks.Filter_for_Tumor_vs_Normal.", run_id, ".tsv")
write.table(x = peaks_filtered_df, file = file2write, sep = "\t", quote = F, row.names = F)
