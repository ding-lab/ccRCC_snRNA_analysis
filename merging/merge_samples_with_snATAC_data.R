# Yige Wu @WashU Jul 2020
## for merging 31 snRNA datasets

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
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
options(future.globals.maxSize = 1000 * 1024^2)

# input dependencies ------------------------------------------------------
## input seurat paths
paths_srat <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20210308.v1/Seurat_Object_Paths.20210308.v1.tsv")
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210305.v1/meta_data.20210305.v1.tsv")

# input the seurat object and store in a list--------------------------------------------------------
paths_srat2process <- paths_srat %>%
  filter(Aliquot %in% idmetadata_df$Aliquot.snRNA[idmetadata_df$snATAC_used])
renal.list <- list()
for (i in 1:nrow(paths_srat2process)) {
  sample_id_tmp <- paths_srat2process$Aliquot[i]
  seurat_obj_path <- paths_srat2process$Path_katmai_seurat_object[i]
  seurat_obj <- readRDS(file = seurat_obj_path)
  renal.list[[sample_id_tmp]] <- seurat_obj
}
length(renal.list)
rm(seurat_obj)

# Run the standard workflow for visualization and clustering ------------
renal.list <- lapply(X = renal.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
cat("Finished building seurat object list!\n")
## integrate without anchor
renal.integrated <- merge(x = renal.list[[1]], y = renal.list[2:length(renal.list)], project = "integrated")
cat("Finished merging seurat object list!\n")
rm(renal.list)
## scale data with all the features
renal.integrated <- SCTransform(renal.integrated, vars.to.regress = c("nCount_RNA","percent.mito"), return.only.var.genes = F)
cat("Finished SCTransform!\n")
## keep it consistant with individual processing pipeline
renal.integrated <- RunPCA(renal.integrated, npcs = 30, verbose = FALSE)
cat("Finished RUNPCA!\n")
renal.integrated <- RunUMAP(renal.integrated, reduction = "pca", dims = 1:30)
cat("Finished RUNUMAP!\n")
renal.integrated <- FindNeighbors(renal.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
cat("Finished FindNeighbors!\n")
renal.integrated <- FindClusters(renal.integrated, resolution = 0.5)
cat("Finished FindClusters!\n")
## save as RDS file
saveRDS(object = renal.integrated, file = paste0(dir_out, "snRNA_Merged.Samples_w_snATAC_data.", run_id, ".RDS"), compress = T)
cat("Finished saving the output!\n")
