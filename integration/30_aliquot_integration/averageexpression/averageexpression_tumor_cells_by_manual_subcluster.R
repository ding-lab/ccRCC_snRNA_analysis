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

## argument: path to the barcode-to-tumorsubcluster table
path_barcode2tumorsubcluster_df <- args[4]
cat(paste0("Path to the barcode-to-tumorsubcluster table: ", path_barcode2tumorsubcluster_df, "\n"))
cat("###########################################\n")

## input the barcode-to-tumorsubcluster table
barcode2tumorsubcluster_df <- fread(input = path_barcode2tumorsubcluster_df, data.table = F)
barcode2tumorsubcluster_df <- as.data.frame(barcode2tumorsubcluster_df)
cat("finish reading the barcode-to-tumorsubcluster table!\n")
cat("###########################################\n")

## input srat
cat(paste0("Start reading the seurat object: ", "\n"))
srat <- readRDS(path_srat)
print("Finish reading the seurat object!\n")
cat("###########################################\n")

## add info to the meta data
metadata_tmp <- barcode2tumorsubcluster_df
metadata_tmp$tumor_exp_subcluster.name <- paste0(metadata_tmp$orig.ident, "_MC", metadata_tmp$manual_cluster_id)
rownames(metadata_tmp) <- metadata_tmp$integrated_barcode
srat@meta.data <- metadata_tmp

## change identification for the cells to be aliquot id
Idents(srat) <- "tumor_exp_subcluster.name"

## run average expression
aliquot.averages <- AverageExpression(srat)
print("Finish running AverageExpression!\n")
cat("###########################################\n")

## write output
write.table(aliquot.averages, file = path_output, quote = F, sep = "\t", row.names = T)
cat("Finished saving the output\n")
cat("###########################################\n")
