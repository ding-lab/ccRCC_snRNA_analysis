# Yige Wu @WashU Sep 2019
## for isolating the non-immune cell clusters and re-do clustering

# library or install.packages-----------------------------------------------------------------
packages = c(
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
  library(package = pkg_name_tmp, character.only = T)
}

# set working directory ---------------------------------------------------
baseD = "/diskmnt/Projects/ccRCC_scratch/"
setwd(baseD)

# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set upstream directories ------------------------------------------------
dir_resources <- "/diskmnt/Projects/ccRCC_scratch/Resources/"
dir_snRNA_processed <- paste0(dir_resources, "snRNA_Processed_Data/")
dir_analysis_results <- paste0(dir_resources, "Analysis_Results/")
dir_recluster_results <- paste0(dir_analysis_results, "recluster/")

# set parameters ----------------------------------------------------------
dir_out_parent <- paste0(dir_recluster_results, "process_nonimmune_cells_on_cluster/")
dir.create(dir_out_parent)
dir_out <- paste0(dir_out_parent, run_id, "/")
dir.create(dir_out)

# Input integrated seurat object -------------------------------
path_int_obj <- paste0(dir_analysis_results, "integrate_seurat_objects/20191021.v1/Renal_Integrated.20191021.v1.RDS")
renal.int.object <- readRDS(path_int_obj)

# get aliquot IDs to process ----------------------------------------------
snRNA_aliquot_ids <- unique(renal.int.object@meta.data$orig.ident)
snRNA_aliquot_ids

# input cell type assignment table ----------------------------------------
dir_cluster2celltype_file <- paste0(dir_snRNA_processed, "Cell_Type_Assignment/")
cluster2celltype_id <- "20191022.v1"
name_cluster2celltype_file <- paste0("ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.", cluster2celltype_id, ".tsv")
path_cluster2celltype_file <- paste0(dir_cluster2celltype_file, name_cluster2celltype_file)
cluster2celltype_tab <- fread(input = path_cluster2celltype_file, data.table = F)

# get nonimmune cluster ids -----------------------------------------------
nonimmune_clusters <- cluster2celltype_tab %>%
  filter(Is_Immune == "No") %>%
  select(Cluster)

nonimmune_clusters <- nonimmune_clusters$Cluster

# input the gene2cell type annotation file --------------------------------
dir_marker_files <- paste0(dir_resources, "Kidney_Markers/")
marker_id <- "20191022.v4"
name_marker_file <- paste0("CC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.", marker_id, ".tsv")
path_marker_file <- paste0(dir_marker_files, name_marker_file) 
#gene2cellType_tab <- fread(input = path_marker_file)

# subset object by non-immune clusters ------------------------------------
renal.int.object <- subset(renal.int.object, idents = nonimmune_clusters)
renal.int.object@meta.data %>% head()

# create list of seurat objects by sample ---------------------------------
renal.list <- list()
for (sample_id_tmp in snRNA_aliquot_ids) {
  seurat_obj_tmp <- subset(renal.int.object, subset =  orig.ident == sample_id_tmp)
  
  DefaultAssay(seurat_obj_tmp) <- "RNA" 
  seurat_obj_tmp@assays$SCT <- NULL
  seurat_obj_tmp@graphs <- list()
  seurat_obj_tmp@neighbors <- list()
  seurat_obj_tmp@reductions <- list()
  for (col_name_tmp in c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.5", "seurat_clusters")) {
    if (col_name_tmp %in% names(seurat_obj_tmp@meta.data)) {
      seurat_obj_tmp@meta.data[[col_name_tmp]] <- NULL
    }
  }
  renal.list[[sample_id_tmp]] <- seurat_obj_tmp
}

#  split the combined object into a list, with each dataset as an element ----------------------------------------
renal.list <- lapply(X = renal.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# create reference anchors  ----------------------------------------------------
renal.anchors <- FindIntegrationAnchors(object.list = renal.list, dims = 1:20)
saveRDS(object = renal.anchors, file = paste0(dir_out, "NonImmune_Anchors.", run_id, ".RDS"))

# Integrate Data -----------------------------------------------------------
renal.int.nonimmune.obj <- IntegrateData(anchorset = renal.anchors, dims = 1:20)


# switch to integrated assay. ---------------------------------------------
DefaultAssay(renal.int.nonimmune.obj) <- "integrated" #only have 3000 features

# Run the standard workflow for visualization and clustering ------------
renal.int.nonimmune.obj <- ScaleData(renal.int.nonimmune.obj, verbose = F) 
renal.int.nonimmune.obj <- RunPCA(renal.int.nonimmune.obj, npcs = 30, verbose = FALSE)
renal.int.nonimmune.obj <- RunUMAP(renal.int.nonimmune.obj, reduction = "pca", dims = 1:20)
renal.int.nonimmune.obj <- FindNeighbors(renal.int.nonimmune.obj, reduction = "pca", dims = 1:20)
renal.int.nonimmune.obj <- FindClusters(renal.int.nonimmune.obj, resolution = 0.5)
saveRDS(object = renal.int.nonimmune.obj, file = paste0(dir_out, "NonImmune_Integrated.", run_id, ".RDS"))

# find DEG ----------------------------------------------------------------
DefaultAssay(renal.int.nonimmune.obj) <- "RNA"

renal.markers <- FindAllMarkers(object = renal.int.nonimmune.obj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
renal.markers %>%
  colnames()
renal.markers <- renal.markers[, c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
write.table(renal.markers, file = paste0(dir_out, "NonImmune.DEGs.Pos.txt"), quote = F, sep = "\t", row.names = F)


