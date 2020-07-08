# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(ggridges)
library(viridis)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input barcode2celltype --------------------------------------------------
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells/20200626.v1/30AliquotIntegration.Barcode2CellType.TumorManualCluster.20200626.v1.tsv")
## get all the downloaded desired output
dir_files <- "./Resources/snRNA_Processed_Data/InferCNV/outputs/Individual.20200305.v1/"
names_files <- list.files(path = dir_files, pattern = "expr.infercnv.dat", recursive = T)

# extract the tumor only CNVs for each sample -----------------------------
for (name_file in names_files) {
  ## get aliquot id
  id_aliquot <- gsub(x = name_file, pattern = "/expr.infercnv.dat", replacement = "")
  ## get barcodes for tumor cells
  barcodes_tumorcells_all <- barcode2celltype_df$individual_barcode[barcode2celltype_df$orig.ident == id_aliquot & barcode2celltype_df$Cell_type.shorter == "Tumor cells"]
  if (length(barcodes_tumorcells_all > 1500)) {
    barcodes_tumorcells <- sample(x = barcodes_tumorcells_all, size = 1500, replace = F)
  }
  ## input cnv file
  path_file <- paste0(dir_files, name_file)
  cnv_df <- fread(data.table = F, input = path_file, sep = "\t")
  ## file columns
  cnv_write_df <- cnv_df[, barcodes_tumorcells]
  rownames(cnv_write_df) <- cnv_df[, "V1"]
  ## write output
  file2write <- paste0(dir_out, id_aliquot, ".tumorcells.expr.infercnv.dat")
  write.table(x = cnv_write_df, file = file2write, quote = F, row.names = T, sep = "\t")
}

