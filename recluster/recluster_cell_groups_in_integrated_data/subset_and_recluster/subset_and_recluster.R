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

## argument 4: path to the cluster to cell type file
path_cluster2celltype_df <- args[4]
cat(paste0("Path to the cluster-to-cell-type file: ", path_cluster2celltype_df, "\n"))
cat("###########################################\n")

## argument 5: name of the cell group to recluster
cellgroup2process <- args[5]
cat(paste0("Name of the cell group to recluster: ", cellgroup2process, "\n"))
cat("###########################################\n")

## input cluster-to-cell-type file
cluster2celltype_df <- fread(input = path_cluster2celltype_df, data.table = F)
cat("Finish reading the cluster-to-cell-type file\n")
cat("###########################################\n")

## get the cluster ids to process
clusters2process <- cluster2celltype_df %>%
  filter(Most_Enriched_Cell_Group == cellgroup2process) %>%
  select(Cluster)
clusters2process <- clusters2process$Cluster
cat(paste0("Clusters to recluster: ", clusters2process, "\n"))
cat("###########################################\n")

## input srat
cat(paste0("Start reading the seurat object: ", "\n\n"))
srat <- readRDS(path_srat)
cat("Finish reading the seurat object!\n\n")
cat("###########################################\n")

# subset object by clusters2process
srat <- subset(srat, idents = clusters2process)
cat("Finish subsetting the seurat object!\n")
cat("###########################################\n")

## get variably expressed genes
cat("Start running FindVariableFeatures\n")
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 3000)
cat("###########################################\n")
## scale data with all the features
cat("Start running ScaleData\n")
srat <- ScaleData(srat, features = rownames(srat@assays$RNA@counts))
cat("###########################################\n")
## keep it consistant with individual processing pipeline
cat("Start running RunPCA\n")
srat <- RunPCA(srat, npcs = 50, verbose = FALSE)
cat("###########################################\n")
cat("Start running RunUMAP\n")
srat <- RunUMAP(srat, reduction = "pca", dims = 1:30)
cat("###########################################\n")
cat("Start running FindNeighbors\n")
srat <- FindNeighbors(srat, reduction = "pca", dims = 1:30, force.recalc = T)
cat("###########################################\n")
cat("Start running FindClusters\n")
srat <- FindClusters(srat, resolution = 0.5)
cat("###########################################\n")
## save output
cat("Start saving the reclustered seurat object\n")
saveRDS(object = srat, file = path_output, compress = T)
cat("###########################################\n")
