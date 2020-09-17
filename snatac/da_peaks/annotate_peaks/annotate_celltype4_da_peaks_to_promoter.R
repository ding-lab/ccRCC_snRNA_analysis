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
peaks2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/unite_cellgroup4_da_peaks/20200915.v1/DA_peaks.chromvar.MergedObj.byCell_group4.20200915.v1.tsv")
## input peak annotation
peak_anno_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Peak_Annotation/peaks_merged_obj_annotation.tsv")

# merge -------------------------------------------------------------------
peak_anno_df <- peak_anno_df %>%
  mutate(region = paste0(chr, ":", start, "-", end))
peaks2celltype_anno_df <- merge(x = peaks2celltype_df, y = peak_anno_df, by.x = c("V1"), by.y = c("region"), all.x = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DA_peaks.chromvar.MergedObj.byCell_group4.Annotated.", run_id, ".tsv")
write.table(x = peaks2celltype_anno_df, file = file2write, sep = "\t", quote = F, row.names = F)
