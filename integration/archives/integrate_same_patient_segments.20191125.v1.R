# Yige Wu @WashU Sep 2019
## for integrating  snRNA datasets belong to the same patient

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# set case ids to be processed --------------------------------------------
case_id <- c("C3N-01200")
case_id <- c("C3N-00733")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191125.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Case %in% case_id) %>%
  filter(FACS == "") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object

# set aliquot ids to be processed -----------------------------------------
snRNA_aliquot_ids <- seurat_summary2process$Aliquot

# Input seurat objects -----------------------------------------------------
renal.list <- list()
for (i in 1:nrow(seurat_summary2process)) {
  sample_id_tmp <- seurat_summary2process$Aliquot[i]
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[i]
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

#  split the combined object into a list, with each dataset as an element ----------------------------------------
renal.list <- lapply(X = renal.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# create reference anchors  ----------------------------------------------------
renal.anchors <- FindIntegrationAnchors(object.list = renal.list, dims = 1:20)
saveRDS(object = renal.anchors, file = paste0(dir_out, case_id, ".Tummor_Segments.Anchors.", run_id, ".RDS"))
rm(renal.list)

# Integrate Data -----------------------------------------------------------
renal.integrated <- IntegrateData(anchorset = renal.anchors, dims = 1:20)
## Error: vector memory exhausted (limit reached?)
rm(renal.anchors)

# switch to integrated assay. ---------------------------------------------
DefaultAssay(renal.integrated) <- "integrated" #only have 3000 features

# Run the standard workflow for visualization and clustering ------------
# The variable features of this assay are automatically set during IntegrateData
renal.integrated <- ScaleData(renal.integrated, verbose = F) 
renal.integrated <- RunPCA(renal.integrated, npcs = 30, verbose = FALSE)
renal.integrated <- RunUMAP(renal.integrated, reduction = "pca", dims = 1:20)
renal.integrated <- FindNeighbors(renal.integrated, reduction = "pca", dims = 1:20)
renal.integrated <- FindClusters(renal.integrated, resolution = 0.5)

saveRDS(object = renal.integrated, file = paste0(dir_out, case_id, ".Tummor_Segments.Integrated.", run_id, ".RDS"))
