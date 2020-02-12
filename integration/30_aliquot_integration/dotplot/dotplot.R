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
print("finish reading the seurat object!\n\n")
cat("###########################################\n")

## make the data frame for the dotplot
## get the genes within the cell type marker table
genes2plot <-  intersect(gene2celltype_df$Gene, srat@assays$RNA@counts@Dimnames[[1]])
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
malignant_markers <- as.vector(gene2celltype_df$Gene[gene2celltype_df$Cell_Type_Group == "Malignant_Nephron_Epithelium"])
genes2plot_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 25) >= 1, "features.plot"])
genes2plot_filtered <- c(genes2plot_filtered, 
                         as.vector(plot_matrix[(plot_matrix$features.plot %in% malignant_markers) & (rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 5) >= 1), "features.plot"]))
genes2plot_filtered <- unique(genes2plot_filtered)

## make the dotplot
cat("###########################################\n")
cat("Dotplot now\n")
p <- DotPlot(object = srat, features = genes2plot_filtered, col.min = 0)
p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
p <- p  + RotatedAxis()
p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey50"),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 15, face = "bold"),
               strip.placement = "outside")
cat("Finished Dotplot\n")
cat("###########################################\n")

## save dotplot
png(file = path_output, width = 4000, height = 1200, res = 150)
print(p)
dev.off()
cat("Finished saving Dotplot\n")
cat("###########################################\n")
