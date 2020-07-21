Yige Wu @WashU Jul 2020

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
path_rds <- "./Resources/Analysis_Results/integration/30_aliquot_integration/subset/subset_tumor_cells/20200309.v1/30_aliquot_integration.20200212.v3.RDS.tumor_cells.20200309.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-to-tumorsubcluster table
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells/20200720.v1/30AliquotIntegration.Barcode2CellType.TumorManualCluster.20200720.v1.tsv", data.table = F)
barcode2tumorsubcluster_df <- as.data.frame(barcode2tumorsubcluster_df)
cat("finish reading the barcode-to-tumorsubcluster table!\n")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
cat("###########################################\n")

# modify srat object ------------------------------------------------------
## get barcodes to process
barcode2tumorsubcluster_df$Aliquot_WU <- mapvalues(x = barcode2tumorsubcluster_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
barcode2tumorsubcluster_df <- barcode2tumorsubcluster_df %>%
  mutate(Name_TumorSubcluster = paste0(Aliquot_WU, "_C", (Id_TumorManualCluster + 1)))
## change meta data
srat@meta.data$Name_TumorSubcluster <- mapvalues(x = rownames(srat@meta.data), from = barcode2tumorsubcluster_df$integrated_barcode, to = as.vector(barcode2tumorsubcluster_filtered_df$Name_TumorSubcluster))
## change ident
Idents(srat) <- "Name_TumorSubcluster"

## run average expression
aliquot.averages <- AverageExpression(srat)
print("Finish running AverageExpression!\n")
cat("###########################################\n")

## write output
file2write <- paste0(dir_out, "AverageExpression_ByManualTumorSubcluster.", run_id, ".tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)
cat("Finished saving the output\n")
cat("###########################################\n")
