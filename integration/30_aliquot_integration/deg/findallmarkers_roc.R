#!/usr/bin/env Rscript

## library
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

# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("At least 3 arguments must be supplied\n1.path to the seurat object\n2.path to the cell type marker table\n3.path to the output\n", call.=FALSE)
}
## argument 1: path to seurat object
path_srat <- args[1]
cat(paste0("Path to the seurat object: ", path_srat, "\n"))
cat("###########################################\n")
## argument 2: path to the cell type marker table
path_gene2celltype_df <- args[2]
cat(paste0("Path to the cell type marker table: ", path_gene2celltype_df, "\n"))
cat("###########################################\n")
## argument 3: path to the output
path_output <- args[3]
cat(paste0("Path to the output file: ", path_output, "\n"))
cat("###########################################\n")

## input cell type marker table
gene2celltype_df <- fread(input = path_gene2celltype_df, data.table = F)
cat("finish reading the cell type marker table!\n")
cat("###########################################\n")

## input srat
cat(paste0("Start reading the seurat object: ", "\n\n"))
srat <- readRDS(path_srat)
cat("Finish reading the seurat object!\n\n")
cat("###########################################\n")

## run FindAllMarkers
cat(paste0("Start runing FindAllMarkers ", "\n\n"))
markers_roc <- FindAllMarkers(object = srat, test.use = "roc", only.pos = T, return.thresh = 0.5)
cat("Finish runing FindAllMarkers\n\n")
cat("###########################################\n")

## filter out markers with small power    
markers_roc <- markers_roc %>%
      filter(power > 0)

## annotate gene to cell type markers
cat(paste0("Start annotating DEGs to cell type markers", "\n\n"))
  markers_roc$Cell_Type_Group <- mapvalues(x = ifelse(markers_roc$gene %in% gene2celltype_df$Gene, markers_roc$gene, NA), from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
  markers_roc$Cell_Type1 <- mapvalues(x = ifelse(markers_roc$gene %in% gene2celltype_df$Gene, markers_roc$gene, NA), from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
  markers_roc$Cell_Type2 <- mapvalues(x = ifelse(markers_roc$gene %in% gene2celltype_df$Gene, markers_roc$gene, NA), from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
  markers_roc$Cell_Type3 <- mapvalues(x = ifelse(markers_roc$gene %in% gene2celltype_df$Gene, markers_roc$gene, NA), from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
cat("Finish annotation\n")
cat("###########################################\n")

## write table
write.table(markers_roc, file = path_output, quote = F, sep = "\t", row.names = F)
cat("Finish writing the output\n")
cat("###########################################\n")
