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
# dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input CNV data
cnv_state_bycell_bygene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/annotate_barcode_cnv/annotate_barcode_with_gene_level_cnv_using_cnv_genes/20220131.v1/CNV_State_By_Gene_By_Barcode.20220131.v1.tsv")
## input the integrated data
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_30_ccRCC_tumorcells/20220404.v1/ccRCC.34samples.Tumorcells.SeuratIntegrated.20220404.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")
## spcify assay
assay_process <- "RNA"
slot_process <- "data"
cat(paste0("Assay: ", assay_process, "\n"))

# pre-process---------------------------------------------------------------
srat@meta.data$barcode.individual = str_split_fixed(string = rownames(srat@meta.data), pattern = "_", n = 2)[,1]
srat@meta.data$barcode.unique = paste0(srat@meta.data$orig.ident, "_", srat@meta.data$barcode.individual)
## process cnv data
cnv_state_bycell_bygene_df = cnv_state_bycell_bygene_df %>%
  dplyr::mutate(barcode.unique = paste0(id_aliquot, '_', barcode_individual))

gene_cna = "MYC"
for (gene_cna in c("MYC", "QKI", "ARID1B")) {
  cnv_tmp_df = cnv_state_bycell_bygene_df[cnv_state_bycell_bygene_df$gene_symbol == gene_cna,]
  srat@meta.data$cna_state = mapvalues(x = srat@meta.data$barcode.unique, from = cnv_tmp_df$barcode.unique, to = as.vector(cnv_tmp_df$cna_state))
  srat@meta.data$cell_group = paste0(srat@meta.data$orig.ident, "_", srat@meta.data$cna_state)
  Idents(srat) <- "cell_group" 
  # run average expression --------------------------------------------------
  aliquot.averages <- AverageExpression(srat, assays = assay_process, slot = slot_process)
  print("Finish running AverageExpression!\n")
  cat("###########################################\n")
  
  # write output ------------------------------------------------------------
  file2write <- paste0(dir_out, "30ccRCCtumorcellreclustered.", "avgexp.", assay_process, ".", slot_process, ".", gene_cna,  ".CNA.", run_id, ".tsv")
  write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)
  cat("Finished saving the output\n")
  cat("###########################################\n")
}


