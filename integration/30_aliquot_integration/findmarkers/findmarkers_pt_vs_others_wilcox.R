#!/usr/bin/env Rscript

## library
packages = c(
  "ggplot2",
  "Seurat",
  "dplyr",
  "plyr",
  "data.table"
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0("No ", pkg_name_tmp, " Installed!"))
  } else {
    print(paste0("", pkg_name_tmp, " Installed!"))
  }
  library(package = pkg_name_tmp, character.only = T, quietly = T)
}
cat("Finish loading libraries!\n")
cat("###########################################\n")

## get the path to the seurat object
args = commandArgs(trailingOnly=TRUE)

## argument: directory to the output
path_output_dir <- args[1]
cat(paste0("Path to the output directory: ", path_output_dir, "\n"))
cat("###########################################\n")

## argument 2: filename for the output file
path_output_filename <- args[2]
cat(paste0("Filename for the output: ", path_output_filename, "\n"))
cat("###########################################\n")
path_output <- paste0(path_output_dir, path_output_filename)

## argument : path to seurat object
path_srat <- args[3]
cat(paste0("Path to the seurat object: ", path_srat, "\n"))
cat("###########################################\n")

## argument : path to the barcode to cell type table
path_barcode2celltype_df <- args[4]
cat(paste0("Path to the barcode to cell type marker table: ", path_barcode2celltype_df, "\n"))
cat("###########################################\n")

## input cell type marker table
barcode2celltype_df <- fread(input = path_barcode2celltype_df, data.table = F)
cat("finish reading the barcode-cell-type table!\n")
cat("###########################################\n")

## input srat
cat(paste0("Start reading the seurat object: ", "\n"))
srat <- readRDS(path_srat)
print("Finish reading the seurat object!\n")
cat("###########################################\n")

## add cell type info into meta data
metadata_tmp <- srat@meta.data
metadata_tmp$integrated_barcode <- rownames(srat@meta.data)
metadata_tmp <- merge(metadata_tmp, barcode2celltype_df, by = c("integrated_barcode"), all.x = T)
rownames(metadata_tmp) <- metadata_tmp$integrated_barcode
srat@meta.data <- metadata_tmp

## change identification for the cells to be cell type group
Idents(srat) <- "Most_Enriched_Cell_Type1"

## run findallmarkers
markers_df <- FindMarkers(object = srat, ident.1 = "Proximal tubule", ident.2 = c("Endothelial cells", "Fibroblasts", "Loop of Henle", "Lymphoid lineage immune cells", "Myeloid lineage immune cells"), test.use = "wilcox", only.pos = T, logfc.threshold = 0)
print("Finish running FindMarkers!\n")
markers_df$gene <- rownames(markers_df)
cat("###########################################\n")

## write output
write.table(markers_df, file = path_output, quote = F, sep = "\t", row.names = F)
cat("Finished saving the output\n")
cat("###########################################\n")
