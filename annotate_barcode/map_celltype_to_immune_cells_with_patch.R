# Yige Wu @WashU Oct 2020
## map the cell type shorter and detailed

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
## input Alla's immune cell type assignment in March formated by Yige
immune_barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_immune_cells/20200626.v1/Barcode2ImmuneCellType.20200626.v1.tsv")
# > nrow(immune_barcode2celltype_df)
# [1] 33004
## inpu alla's new immune cell types
immune_barcode2celltype_patch_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Integration_AllClusters/ccRCC_metadata_AK_v20200818.txt")

# map by integrated barcodes ----------------------------------------------
## make new cell type shorter
table(immune_barcode2celltype_patch_df$Cell_type.shorter)
celltype.shorter_new_vec <- mapvalues(x = immune_barcode2celltype_df$integrated_barcode, from = immune_barcode2celltype_patch_df$V1, to = as.vector(immune_barcode2celltype_patch_df$Cell_type.shorter))
table(celltype.shorter_new_vec)
## make new cell type detailed
celltype.detailed_new_vec <- mapvalues(x = immune_barcode2celltype_df$integrated_barcode, from = immune_barcode2celltype_patch_df$V1, to = as.vector(immune_barcode2celltype_patch_df$Cell_type.detailed))
table(celltype.detailed_new_vec)

# add in the new cell type shorter and detailed ---------------------------
immune_barcode2celltype_new_df <- immune_barcode2celltype_df 
immune_barcode2celltype_new_df$Cell_type.shorter <- celltype.shorter_new_vec
immune_barcode2celltype_new_df$Cell_type.detailed <- celltype.detailed_new_vec

# write outputs -----------------------------------------------------------
nrow(immune_barcode2celltype_new_df)
# [1] 33004
file2write <- paste0(dir_out, "Barcode2ImmuneCellType.", run_id, ".tsv")
write.table(x = immune_barcode2celltype_new_df, file = file2write, sep = '\t', quote = F, row.names = F)


