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

## argument 1: path to the output directory
path_output_dir <- args[1]
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

## argument: the barcode to cell type table
path_barcode2celltype_df <- args[4]
cat(paste0("Path to the barcode-cell-type info: ",path_barcode2celltype_df, "\n"))
cat("###########################################\n")

## argument : path to the cell type marker table
path_gene2group_df <- args[5]
cat(paste0("Path to the gene list: ", path_gene2group_df, "\n"))
cat("###########################################\n")

## argument: the column name for the gene symbol
genesymbol_col_name <- args[6]
cat(paste0("The column name for the gene symbol column:", genesymbol_col_name, "\n"))
cat("###########################################\n")

## argument: the column name for the gene group
genegroup_col_name <- args[7]
cat(paste0("The column name for the gene group:", genegroup_col_name, "\n"))
cat("###########################################\n")

## argument : minimal percentage of cell in any clustering expressing the genes to show
min.exp.pct <- args[8]
cat(paste0("Minimal percentage of cell in any clustering expressing the genes to show: ", min.exp.pct, "%\n"))
cat("###########################################\n")

## input cell type marker table
gene2group_df <- fread(input = path_gene2group_df, data.table = F)
cat("Finish reading the cell type marker table!\n")
cat("###########################################\n")

## input srat
cat(paste0("Start reading the seurat object: ", "\n\n"))
srat <- readRDS(path_srat)
print("Finish reading the seurat object!\n\n")
cat("###########################################\n")

## make the data frame for the dotplot
## get the genes within the cell type marker table
genes2plot <-  intersect(gene2group_df[, genesymbol_col_name], srat@assays$RNA@counts@Dimnames[[1]])
genes2plot <- unique(genes2plot)
## get the pct expressed for each gene in each cluster
p <- DotPlot(object = srat, features = genes2plot, col.min = 0)
plot_data <- p$data
## transform the dataframe to matrix to better filter out genes with too low expressin
library(data.table)
plot_matrix <- dcast(data = plot_data, formula = features.plot ~ id, value.var = "pct.exp")
plot_matrix %>% head()
## filter for genes that are expressed in >25% of one cluster at least
## replot with the filtered genes plus malignant cell marker genes
genes2plot_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > min.exp.pct) >= 1, "features.plot"])
genes2plot_filtered <- unique(genes2plot_filtered)

## make the dotplot
cat("###########################################\n")
cat("Dotplot now\n")
p <- DotPlot(object = srat, features = genes2plot_filtered, col.min = 0, group.by = "orig.ident")
cat("Finished Dotplot\n")
cat("###########################################\n")

## save dotplot
saveRDS(object = p, file = path_output, compress = T)
cat("Finished saving the output\n")
cat("###########################################\n")
