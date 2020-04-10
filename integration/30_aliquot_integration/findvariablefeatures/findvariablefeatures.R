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

## input srat
cat(paste0("Start reading the seurat object: ", "\n"))
srat <- readRDS(path_srat)
print("Finish reading the seurat object!\n")
cat("###########################################\n")

## get variably expressed genes
cat("Start running FindVariableFeatures\n")
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 3000)
cat("###########################################\n")

## write output
var_features_df <- data.frame(gene = srat@assays$RNA@var.features)
write.table(var_features_df, file = path_output, quote = F, sep = "\t", row.names = F)
cat("Finished saving the output\n")
cat("###########################################\n")
