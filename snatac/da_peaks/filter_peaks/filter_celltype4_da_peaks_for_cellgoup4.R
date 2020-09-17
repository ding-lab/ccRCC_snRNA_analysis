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
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/filter_markers/filter_markers_wilcox_bygroup/20200916.v1/findallmarkers_wilcox_bycellgroup.pos.logfcthreshold0.1.minpct0.1.mindiffpct0.1.Top50avg_logFC.tsv")

# filter DEGs -------------------------------------------------------------
degs_filtered_df <- degs_df %>%
  filter(p_val_adj < 0.05) %>%
  filter(cluster != "Unknown")

# filter by promoter region and DEGs --------------------------------------
peaks_filtered_df <- peaks2celltype_df %>%
  filter(avg_logFC > 0) %>%
  filter(Cell_type.filename != "Unknown") %>%
  filter(annotation == "Promoter")
## filter DEGs
peaks_degs_df <- merge(x = peaks_filtered_df, y = degs_filtered_df, by.x = c("SYMBOL", "Cell_type.filename"), by.y = c("gene", "cluster"), suffixes = c(".peak", ".deg"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cell_Type_Specific_DEGs_in_DA_peaks.", run_id, ".tsv")
write.table(x = peaks_degs_df, file = file2write, sep = "\t", quote = F, row.names = F)
