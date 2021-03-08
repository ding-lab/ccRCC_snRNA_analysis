#!/usr/bin/env Rscript

## library
cat("###########################################\n")
packages = c(
  "ggplot2",
  "Seurat",
  "dplyr",
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

## argument 1: path to the output directory
path_output_dir <- args[1]
dir.create(path_output_dir)
cat(paste0("Path to the output directory: ", path_output_dir, "\n"))
cat("###########################################\n")

## argument 2: filename for the output file
path_output_filename <- args[2]
cat(paste0("Filename for the output: ", path_output_filename, "\n"))
cat("###########################################\n")
path_output <- paste0(path_output_dir, path_output_filename)

## argument 3: path to seurat object
path_srat <- args[3]
cat(paste0("Path to the seurat object: ", path_srat, "\n"))
cat("###########################################\n")

## input srat
cat(paste0("Start reading the seurat object: ", "\n"))
srat <- readRDS(path_srat)
cat("Finish reading the seurat object!\n")
cat("###########################################\n")

## fetch data
cat("Start fetching the data\n")
umap_data <- Seurat::FetchData(object = srat, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
umap_data$barcode <- rownames(umap_data)
cat("###########################################\n")
cat("Start writing the output\n")
write.table(x = umap_data, file = path_output, quote = F, row.names = F, sep = "\t")
cat("###########################################\n")
