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
  stop("At least 3 arguments must be supplied\n1.path to the seurat object\n2.path to the cell type marker table\n3.path to the output directory\n", call.=FALSE)
}
## argument 1: path to seurat object
path_srat <- args[1]
cat(paste0("Path to the seurat object: ", path_srat, "\n"))
cat("###########################################\n")
## argument 2: path to the cell type marker table
path_gene2celltype_df <- args[2]
cat(paste0("Path to the cell type marker table: ", path_gene2celltype_df, "\n"))
cat("###########################################\n")
## argument 3: directory to the output
path_output_dir <- args[3]
cat(paste0("Path to the output directory: ", path_output_dir, "\n"))
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

## create directories for cell groups
cellgroups <- gene2celltype_df %>%
  select(Cell_Type_Group) %>%
  unique()
for (i in 1:nrow(cellgroups)) {
  dir2create <- paste0(c(path_output_dir, cellgroups[i,]), collapse = "/")
  dir.create(dir2create)
}

## create directories for cell type 1
cellgroups <- gene2celltype_df %>%
  select(Cell_Type_Group, Cell_Type1) %>%
  unique()
for (i in 1:nrow(cellgroups)) {
  dir2create <- paste0(c(path_output_dir, unlist(cellgroups[i,])), collapse = "/")
  dir.create(dir2create)
}

## create directories for cell type 2
cellgroups <- gene2celltype_df %>%
  select(Cell_Type_Group, Cell_Type1, Cell_Type2) %>%
  unique()
for (i in 1:nrow(cellgroups)) {
  dir2create <- paste0(c(path_output_dir, unlist(cellgroups[i,])), collapse = "/")
  dir.create(dir2create)
}

## create directories for cell type 3
cellgroups <- gene2celltype_df %>%
  select(Cell_Type_Group, Cell_Type1, Cell_Type2, Cell_Type3) %>%
  unique()
for (i in 1:nrow(cellgroups)) {
  dir2create <- paste0(c(path_output_dir, unlist(cellgroups[i,])), collapse = "/")
  dir.create(dir2create)
}

## create directories for cell type 4
cellgroups <- gene2celltype_df %>%
  select(Cell_Type_Group, Cell_Type1, Cell_Type2, Cell_Type3, Cell_Type4) %>%
  unique()
for (i in 1:nrow(cellgroups)) {
  dir2create <- paste0(c(path_output_dir, unlist(cellgroups[i,])), collapse = "/")
  dir.create(dir2create)
}

## make the dotplot
cat("###########################################\n")
for (i in 1:nrow(gene2celltype_df)) {
  gene_tmp <- gene2celltype_df[i,"Gene"]
    ## get path to the output diretory
    path_output_dir_tmp <- paste0(c(path_output_dir, unlist(gene2celltype_df[i, c("Cell_Type_Group", "Cell_Type1", "Cell_Type2", "Cell_Type3", "Cell_Type4")])), collapse = "/")
    path_output <- paste0(path_output_dir_tmp, gene_tmp, ".30_aliquot_integration.png")
    if (!file.exists(path_output) & (gene_tmp %in% srat@assays$RNA@counts@Dimnames[[1]]) {
    ## make featureplot
      p <- FeaturePlot(object = srat, features = gene_tmp, 
        min.cutoff = "q10", max.cutoff = "q90", sort.cell = TRUE,
                       cols = c("grey", "red"), reduction = "umap", label = T)
    
    ## save featureplot
    png(file = path_output, width = 4000, height = 4000, res = 150)
    print(p)
    dev.off()
}
}
