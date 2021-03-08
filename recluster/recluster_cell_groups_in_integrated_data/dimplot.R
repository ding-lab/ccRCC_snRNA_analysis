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

## argument 4: format of the figure
format_figure <- args[4]
cat(paste0("Format of the figure: ", format_figure, "\n"))
cat("###########################################\n")

## argument 5: width of the figure
width <- args[5]
cat(paste0("Width of the figure: ", width, "\n"))
cat("###########################################\n")

## argument 6: height of the figure
height <- args[6]
cat(paste0("Height of the figure: ", height, "\n"))
cat("###########################################\n")

## input srat
cat(paste0("Start reading the seurat object: ", "\n"))
srat <- readRDS(path_srat)
cat("Finish reading the seurat object!\n")
cat("###########################################\n")

## make the dimplot
cat("Start making the dimplot\n")
p <- DimPlot(object = srat, label = TRUE, label.size = 10) + NoLegend()
cat("###########################################\n")

## save output
cat("Start saving the object\n")
if (format_figure == "PNG") {
    png(file = path_output, width = as.numeric(width), height = as.numeric(height), res = 150)
    print(p)
    dev.off()
}
if (format_figure == "PDF") {
pdf(file = path_output, width = height, height = height)
print(p)
dev.off()
}
cat("###########################################\n")
