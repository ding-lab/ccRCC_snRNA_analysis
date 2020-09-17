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
path_rds <- "./Resources/Analysis_Results/integration/31_aliquot_integration/31_aliquot_integration_without_anchoring/20200727.v1/31_aliquot_integration_without_anchoring.20200727.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20200910.v2/31Aliquot.Barcode2CellType.20200910.v2.tsv", data.table = F)
barcode2celltype_df <- as.data.frame(barcode2celltype_df)
cat("finish reading the barcode-to-cell type table!\n")
## spcify assay
assay_process <- "SCT"
cat(paste0("Assay: ", assay_process, "\n"))

# set ident ---------------------------------------------------------------
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_bycellgroup_byaliquot = paste0(orig.ident, "_", Cell_group3)) %>%
  mutate(id_cell = paste0(orig.ident, "_", individual_barcode))
srat@meta.data$individual_barcode <- str_split_fixed(string = rownames(srat@meta.data), pattern = "_", n = 2)[,1]
srat@meta.data$id_cell <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$individual_barcode)
srat@meta.data$id_bycellgroup_byaliquot <- mapvalues(x = srat@meta.data$id_cell, from = barcode2celltype_df$id_cell, to = as.vector(barcode2celltype_df$id_bycellgroup_byaliquot))
Idents(srat) <- "id_bycellgroup_byaliquot" 

# run average expression --------------------------------------------------
aliquot.averages <- AverageExpression(srat, assays = assay_process, slot = "scale.data")
print("Finish running AverageExpression!\n")
cat("###########################################\n")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "avgexp.", assay_process, ".bycellgroup3.byaliquot.", "31_aliquot_integration.", run_id, ".tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)
cat("Finished saving the output\n")
cat("###########################################\n")


