# Yige Wu @WashU Aug 2020
## for integrating  snRNA datasets belong to the same patient

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
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
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

# set case ids to be processed --------------------------------------------
ids_case_process <- c("C3L-00416", "C3N-00242")

# input dependencies ------------------------------------------------
## input srat paths
paths_srat_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20200717.v1/Seurat_Object_Paths.20200717.v1.tsv")
## input barcode-to-celltype info
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells/20200720.v1/30AliquotIntegration.Barcode2CellType.TumorManualCluster.20200720.v1.tsv", data.table = F)

# set aliquot ids to be processed -----------------------------------------
paths_srat_df <- paths_srat_df %>%
  filter(Case %in% ids_case_process)

# Input seurat objects -----------------------------------------------------
renal.list <- list()
for (i in 1:nrow(paths_srat_df)) {
  ## input seurat object
  sample_id_tmp <- paths_srat_df$Aliquot[i]
  seurat_obj_path <- paths_srat_df$Path_katmai_seurat_object[i]
  seurat_obj <- readRDS(file = seurat_obj_path)
  print(dim(seurat_obj))
  
  ## subset to tumor cells
  ### get tumor cell barcode
  barcodes_tumorcells <- barcode2celltype_df$individual_barcode[barcode2celltype_df$orig.ident == sample_id_tmp & barcode2celltype_df$Cell_type.shorter == "Tumor cells"]
  seurat_obj <- subset(seurat_obj, cells = barcodes_tumorcells)
  print(dim(seurat_obj))
  
  ## clean
  seurat_obj$orig.ident  <- sample_id_tmp
  DefaultAssay(seurat_obj) <- "RNA" 
  seurat_obj@assays$SCT <- NULL
  seurat_obj@graphs <- list()
  seurat_obj@neighbors <- list()
  seurat_obj@reductions <- list()
  for (col_name_tmp in c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.5", "seurat_clusters")) {
    if (col_name_tmp %in% names(seurat_obj@meta.data)) {
      seurat_obj@meta.data[[col_name_tmp]] <- NULL
    }
  }
  renal.list[[sample_id_tmp]] <- seurat_obj
  
}
length(renal.list)
rm(seurat_obj)

#  split the combined object into a list, with each dataset as an element ----------------------------------------
renal.list <- lapply(X = renal.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = num_features)
})

# Merge Data -----------------------------------------------------------
renal.merged <- merge(x = renal.list[[1]], y = renal.list[2:length(renal.list)], project = "Merged")

# Run the standard workflow for visualization and clustering ------------
# The variable features of this assay are automatically set during IntegrateData
renal.merged <- SCTransform(renal.merged, vars.to.regress = c("nCount_RNA","percent.mito"), return.only.var.genes = F)
renal.merged <- RunPCA(renal.merged, npcs = num_pc, verbose = FALSE)
renal.merged <- RunUMAP(renal.merged, reduction = "pca", dims = 1:num_pc)
renal.merged <- FindNeighbors(renal.merged, reduction = "pca", dims = 1:num_pc)
renal.merged <- FindClusters(renal.merged, resolution = resolution_findclusters)
saveRDS(object = renal.merged, file = paste0(dir_out, paste0(ids_case_process, collapse = "_"), ".Tumor_Segments.Merged.", run_id, ".RDS"))

