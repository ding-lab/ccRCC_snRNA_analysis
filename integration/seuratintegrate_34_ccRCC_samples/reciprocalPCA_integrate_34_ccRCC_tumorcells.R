# Yige Wu @WashU Mar 2021
## https://satijalab.org/seurat/archive/v3.0/integration.html
## also used references

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
# dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "future",
  "future.apply"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
# set up future for parallelization
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
## input seurat paths
paths_srat <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20210729.v1/Seurat_Object_Paths.20210729.v1.tsv")
## input the cell to cell type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
## input doublet prediction
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)

# preprocess ------------------------------------------------------------------
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_sample_barcode = paste0(orig.ident, "_", individual_barcode))
srat@meta.data$individual_barcode <- str_split_fixed(string = rownames(srat@meta.data), pattern = "_", n = 2)[,1]
srat@meta.data$id_sample_barcode <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$individual_barcode)
srat@meta.data$Cell_type.shorter <- mapvalues(x = srat@meta.data$id_sample_barcode, from = barcode2celltype_df$id_sample_barcode, to = as.vector(barcode2celltype_df$Cell_type.shorter))
barcodes_keep <- rownames(srat@meta.data)[srat@meta.data$Cell_type.shorter %in% c("Tumor cells", "EMT tumor cells")]

# process -----------------------------------------------------------------
path_anchor_file <- paste0(dir_out, "anchor" , ".RDS")

if (!file.exists(path_anchor_file)) {
  paths_srat2process <- paths_srat %>%
    filter(Case != "C3L-00359")
  srat_list <- list()
  for (i in 1:nrow(paths_srat2process)) {
    sample_id_tmp <- paths_srat2process$Aliquot[i]
    seurat_obj_path <- paths_srat2process$Path_katmai_seurat_object[i]
    seurat_obj <- readRDS(file = seurat_obj_path)
    ## subset the barcodes
    barcode2scrublet_tmp_df <- barcode2scrublet_df %>%
      filter(Aliquot == sample_id_tmp) %>%
      filter(Barcode %in% rownames(seurat_obj@meta.data)) %>%
      filter(!predicted_doublet)
    barcodes_keep <- barcode2scrublet_tmp_df$Barcode; length(barcodes_keep)
    barcodes_tumorcells <- barcode2celltype_df$individual_barcode[barcode2celltype_df$Cell_type.shorter %in% c("Tumor cells", "EMT tumor cells")]
    barcodes_keep <- barcodes_keep[(barcodes_keep %in% barcodes_tumorcells)]; length(barcodes_keep)
    ###
    print("subsetting")
    seurat_sub_obj <- subset(x = seurat_obj, cells = barcodes_keep)
    print(dim(seurat_sub_obj))
    srat_list[[i]] <- seurat_sub_obj
  }
  length(srat_list)
  rm(seurat_obj)
  cat("Finished making seurat object list!\n")
  
  # Run the standard workflow for visualization and clustering ------------
  srat_list <- future_lapply(X = srat_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })
  cat("Finished NormalizeData and FindVariableFeatures for the list!\n")
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = srat_list)
  cat("Finished SelectIntegrationFeatures!\n")
  
  srat_list <- future_lapply(X = srat_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  ## chose one male and one female as references to scale up the integration
  ### C3N-00437	 is female
  ### C3L-00088 is male
  srat_anchors <- FindIntegrationAnchors(object.list = srat_list, anchor.features = features, reduction = "rpca", dims = 1:30, reference = c(3, 5))
  cat("Finished FindIntegrationAnchors!\n")
  
  saveRDS(object = srat_anchors, file = path_anchor_file, compress = T)
  cat("Finished saveRDS srat_anchors!\n")
}  else {
  srat_anchors <- readRDS(file = path_anchor_file)
}

# this command creates an 'integrated' data assay
srat_integrated <- IntegrateData(anchorset = srat_anchors, dims = 1:30)
cat("Finished IntegrateData!\n")
rm(srat_list)
rm(srat_anchors)
cat("Finished deleting list and anchor!\n")

## keep it consistant with individual processing pipeline
DefaultAssay(srat_integrated)
# DefaultAssay(pancreas.integrated) <- "integrated"
srat_integrated <- ScaleData(srat_integrated, verbose = F)
cat("Finished ScaleData!\n")
# srat_integrated <- RunPCA(srat_integrated, npcs = 30, verbose = FALSE)
srat_integrated <- RunPCA(srat_integrated, npcs = 45, verbose = FALSE)
cat("Finished RUNPCA!\n")
srat_integrated <- RunUMAP(srat_integrated, reduction = "pca", dims = 1:30)
cat("Finished RUNUMAP!\n")
srat_integrated <- FindNeighbors(srat_integrated, reduction = "pca", dims = 1:30, force.recalc = T)
cat("Finished FindNeighbors!\n")
srat_integrated <- FindClusters(srat_integrated, resolution = 0.5)
cat("Finished FindClusters!\n")
## save as RDS file
## scale data with all the features
path_final_file <- paste0(dir_out, "ccRCC.34samples.Tumorcells.SeuratIntegrated.", run_id, ".RDS")
saveRDS(object = srat_integrated, file = path_final_file, compress = T)
cat("Finished saving the srat_integrated!\n")

# fetch data --------------------------------------------------------------
umap_data <- Seurat::FetchData(object = srat_integrated, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
umap_data$barcode <- rownames(umap_data)

file2write <- paste0(dir_out, "ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.", run_id, ".tsv")
write.table(x = umap_data, file = file2write, quote = F, sep = "\t", row.names = F)
print("Finish writing the output!\n")

# plot --------------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.", run_id, ".png")
png(filename = file2write, width = 1200, height = 1100, res = 150)
DimPlot(srat_integrated,reduction = "umap")
dev.off()

