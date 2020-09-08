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
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_immune_cellgroups/20200908.v1/31Aliquot.Barcode2CellType.Immune_Grouped.20200908.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")
## spcify assay
assay_process <- "SCT"
cat(paste0("Assay: ", assay_process, "\n"))

# set ident ---------------------------------------------------------------
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_cell = paste0(orig.ident, "_", individual_barcode))
srat@meta.data$individual_barcode <- str_split_fixed(string = rownames(srat@meta.data), pattern = "_", n = 2)[,1]
## check if the individual_barcode is mapped right
srat@meta.data %>% head()
## add cell id to the seurat meta data
srat@meta.data$id_cell <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$individual_barcode)
srat@meta.data$Cell_group.immune <- mapvalues(x = srat@meta.data$id_cell, from = barcode2celltype_df$id_cell, to = as.vector(barcode2celltype_df$Cell_group.immune))
unique(srat@meta.data$Cell_group.immune)
Idents(srat) <- "Cell_group.immune" 
# run average expression --------------------------------------------------
aliquot.averages <- AverageExpression(srat, assays = assay_process, slot = "scale.data")
print("Finish running AverageExpression!\n")
cat("###########################################\n")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "averageexpression.", assay_process, ".byCell_group.immune.", "31_aliquot_integration.", run_id, ".tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)
cat("Finished saving the output\n")
cat("###########################################\n")


