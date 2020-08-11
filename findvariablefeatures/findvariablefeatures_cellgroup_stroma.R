# Yige Wu @WashU Aug 2020

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
library(tibble)
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
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200811.v1/31Aliquot.Barcode2CellType.20200811.v1.tsv", data.table = F)
barcode2celltype_df <- as.data.frame(barcode2celltype_df)
cat("finish reading the barcode-to-cell type table!\n")
## specify the cell group
cellgroup2process <- "Stroma"

# subset seurat object ----------------------------------------------------
BC <- srat@meta.data %>% rownames
#### get original barcode
srat@meta.data$original_barcode <- BC %>% strsplit("_") %>% lapply("[[",1) %>% unlist
srat@meta.data <- merge(srat@meta.data,barcode2celltype_df,
                            by.x=c("orig.ident","original_barcode"), by.y=c("orig.ident","individual_barcode"))
srat@meta.data <- column_to_rownames(srat@meta.data,"rowname")[BC,]
Idents(srat) <- "Cell_group.shorter"
## subset data
cat("Start subsetting\n")
srat.new <- subset(srat, idents = "Stroma")
cat(dim(srat.new))
rm(srat)
## set default assay
DefaultAssay(object = srat.new) <- "RNA"
cat("###########################################\n")

# find variable genes -----------------------------------------------------
cat("Start running FindVariableFeatures\n")
srat.new <- FindVariableFeatures(srat.new, selection.method = "vst", nfeatures = 3000, verbose = T)
cat("###########################################\n")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "findvariablefeatures.", "cellgroup.", cellgroup2process, ".", run_id, ".tsv")
var_features_df <- data.frame(gene = srat.new@assays$RNA@var.features)
write.table(var_features_df, file = file2write, quote = F, sep = "\t", row.names = F)
cat("Finished saving the output\n")
cat("###########################################\n")


