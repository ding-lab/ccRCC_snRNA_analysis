# Yige Wu @WashU Jul 2020
## for integrating  snRNA datasets with normal epithelial cells

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
paths_srat <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20200717.v1/Seurat_Object_Paths.20200717.v1.tsv")

# set case ids to be processed --------------------------------------------
ids_case2process <- c("C3L-00088", "C3N-01200", "C3N-01213")
## get the srat object to process
paths_srat2process <- paths_srat %>%
  filter(Case %in% ids_case2process)

# make a list of seurat objects -----------------------------------------------------
renal.list <- list()
for (i in 1:nrow(paths_srat2process)) {
  sample_id_tmp <- paths_srat2process$Aliquot[i]
  seurat_obj_path <- paths_srat2process$Path_box_seurat_object[i]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
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

# Merge Data -----------------------------------------------------------
renal.merged <- merge(renal.list[[1]], renal.list[2:length(renal.list)], project = "Merged")

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, case_id, ".Tumor_Segments.Merged.", run_id, ".RDS")
saveRDS(object = renal.merged, file = file2write, compress = T)

