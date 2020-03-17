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

## argument : path to the cell type marker table
path_gene2celltype_df <- args[4]
cat(paste0("Path to the cell type marker table: ", path_gene2celltype_df, "\n"))
cat("###########################################\n")

## input cell type marker table
gene2celltype_df <- fread(input = path_gene2celltype_df, data.table = F)
cat("finish reading the cell type marker table!\n")
cat("###########################################\n")

## input srat
cat(paste0("Start reading the seurat object: ", "\n"))
srat <- readRDS(path_srat)
print("Finish reading the seurat object!\n")
cat("###########################################\n")

## run findallmarkers
markers_roc <- FindAllMarkers(object = srat, test.use = "roc", only.pos = T, return.thresh = 0.5)
print("Finish running FindAllMarkers!\n")
cat("###########################################\n")

## filter by clster distinguishing power
markers_roc <- markers_roc %>%
      filter(power > 0)

## annotate genes to cell types
markers_roc$Cell_Type_Group <- mapvalues(x = ifelse(markers_roc$gene %in% gene2celltype_df$Gene, markers_roc$gene, NA), from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
markers_roc$Cell_Type1 <- mapvalues(x = ifelse(markers_roc$gene %in% gene2celltype_df$Gene, markers_roc$gene, NA), from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
markers_roc$Cell_Type2 <- mapvalues(x = ifelse(markers_roc$gene %in% gene2celltype_df$Gene, markers_roc$gene, NA), from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
markers_roc$Cell_Type3 <- mapvalues(x = ifelse(markers_roc$gene %in% gene2celltype_df$Gene, markers_roc$gene, NA), from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)

## write output
write.table(markers_roc, file = path_output, quote = F, sep = "\t", row.names = F)
cat("Finished saving the output\n")
cat("###########################################\n")
