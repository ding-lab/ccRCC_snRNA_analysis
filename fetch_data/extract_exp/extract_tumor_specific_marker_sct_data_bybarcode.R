# Yige Wu @WashU Sep 2020
# Yige Wu @WashU Jun 2021
## use RNA assay according to https://github.com/satijalab/seurat/issues/2646

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 5 * 1024^3) # for 5 Gb RAM

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Data_Freezes/V2/snRNA/All_Cells_Merged/33_aliquot_merged_without_anchoring.20210428.v2.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Data_Freezes/V2/snRNA/Cell_Type_Assignment/33Aliquot.Barcode2CellType.20210423.v1.tsv", data.table = F)
## input the markers
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210702.v1/ccRCC_markers.Surface.20210702.v1.tsv")

# preprocess the Seurat object meta data---------------------------------------------
BC <- srat@meta.data %>% rownames
## get original barcode
srat@meta.data$original_barcode <- BC %>% strsplit("_") %>% lapply("[[",1) %>% unlist
head(srat@meta.data$original_barcode)
cat("finish adding the simple barcode!\n")
## make combined id for the seurat meta data
srat@meta.data$id_aliquot_barcode <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$original_barcode)
# head(srat@meta.data$id_aliquot_barcode)
ids_aliquot_barcode_currentall <- srat@meta.data$id_aliquot_barcode
cat("finish adding unique id for each barcode in the seurat object!\n")

# preprocess barcode info -------------------------------------------------
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode))
ids_aliquot_barcode_process <- barcode2celltype_df$id_aliquot_barcode[barcode2celltype_df$Cell_group5 == "Tumor cells"]
barcodes_process <- BC[ids_aliquot_barcode_currentall %in% ids_aliquot_barcode_process]
ids_aliquot_barcode_process <- ids_aliquot_barcode_currentall[ids_aliquot_barcode_currentall %in% ids_aliquot_barcode_process]
barcode_map_df <- data.frame(barcode_merged = barcodes_process, aliquot_barcode = ids_aliquot_barcode_process)

# extract tumor-cell expression of the given genes ------------------------
genes_process <- genes_process_df$Gene
DefaultAssay(srat) <- "SCT"
exp_df <- FetchData(object = srat, vars = genes_process, slot = "data")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "tumor_specific_marker_expression_bybarcode.RDS")
saveRDS(object = exp_df, file = file2write, compress = T)
file2write <- paste0(dir_out, "barcde_mapping.tsv")
write.table(x = barcode_map_df, file = file2write, sep = "\t", row.names = F, quote = F)


