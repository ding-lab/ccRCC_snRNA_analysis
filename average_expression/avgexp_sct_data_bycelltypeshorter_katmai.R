# Yige Wu @WashU Aug 2020
## finding differentially expressed gene for each cell type using integrared object

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

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Data_Freezes/V2/snRNA/All_Cells_Merged/33_aliquot_merged_without_anchoring.20210428.v2.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Data_Freezes/V2/snRNA/Cell_Type_Assignment/33Aliquot.Barcode2CellType.20210423.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")
## spcify assay
assay_process <- "SCT"
slot_process <- "data"
cat(paste0("Assay: ", assay_process, "\n"))

# set ident ---------------------------------------------------------------
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_cell = paste0(orig.ident, "_", individual_barcode))
srat@meta.data$individual_barcode <- str_split_fixed(string = rownames(srat@meta.data), pattern = "_", n = 2)[,1]
## check if the individual_barcode is mapped right
srat@meta.data %>% head()
## add cell id to the seurat meta data
srat@meta.data$id_cell <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$individual_barcode)
srat@meta.data$Cell_group <- mapvalues(x = srat@meta.data$id_cell, from = barcode2celltype_df$id_cell, to = as.vector(barcode2celltype_df$Cell_type.shorter))
unique(srat@meta.data$Cell_group)
Idents(srat) <- "Cell_group" 

# run average expression --------------------------------------------------
aliquot.averages <- AverageExpression(srat, assays = assay_process, slot = slot_process)
print("Finish running AverageExpression!\n")
cat("###########################################\n")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "33_aliquot_merged.", "avgexp.", assay_process, ".", slot_process, ".", "byCell_type.shorter.", run_id, ".tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)
cat("Finished saving the output\n")
cat("###########################################\n")


