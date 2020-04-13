# Yige Wu @WashU Apr 2020
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
path_rds <- "./Resources/Analysis_Results/integration/30_aliquot_integration/docker_run_integration/20200212.v3/30_aliquot_integration.20200212.v3.RDS"
srat <- readRDS(file = path_rds)
## input the barcode-cell-type table
path_barcode2celltype <- "./Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200410.v1/30_aliquot_integration.barcode2celltype.20200410.v1.tsv"
barcode2celltype_df <- fread(input = path_barcode2celltype, data.table = F)
## input id meta data
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## set tumor aliquot id
id_aliquot <- "CPT0075130004"
id_case <- id_metadata_df$Case[id_metadata_df$Aliquot.snRNA == id_aliquot]
id_case

# set ident ---------------------------------------------------------------
barcode2celltype_df$id_case <- mapvalues(x = barcode2celltype_df$orig.ident, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))
## get the barcode name for tumor cells and normal proximal tubule
integratedbarcodes_group1 <- barcode2celltype_df$integrated_barcode[barcode2celltype_df$orig.ident == id_aliquot & barcode2celltype_df$Cell_type.detailed == "Tumor cells"]
integratedbarcodes_group2 <- barcode2celltype_df$integrated_barcode[barcode2celltype_df$id_case == id_case & barcode2celltype_df$Cell_type.detailed == "Proximal tubule"]
## modify meta data
srat@meta.data <- barcode2celltype_df
rownames(srat@meta.data) <- barcode2celltype_df$integrated_barcode
# run findmarkers ------------------------------------------------------
marker_roc_df <- FindMarkers(object = srat, test.use = "roc", return.thresh = 0.5, verbose = T, cells.1 = integratedbarcodes_group1, cells.2 = integratedbarcodes_group2)
marker_roc_df$row_name <- rownames(marker_roc_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "findmarkers_roc.", id_aliquot, "tumor_vs_normalpt.",  run_id, ".tsv")
write.table(x = marker_roc_df, file = file2write, sep = "\t", quote = F, row.names = F)


