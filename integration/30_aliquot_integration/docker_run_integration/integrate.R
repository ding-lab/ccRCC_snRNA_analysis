# Yige Wu @WashU Sep 2019
## for integrating two snRNA datasets for sample CPT0086820004 and CPT0075130004 (from cellranger output with premrna reference)

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
dir_base = "/diskmnt/Projects/ccRCC_scratch/"
setwd(dir_base)

# set run id --------------------------------------------------------------
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
print(run_id)


# set upstream directories ------------------------------------------------
dir_resources <- paste0(dir_base, "Resources/")
dir_analysis_results <- paste0(dir_resources, "Analysis_Results/")
dir_snRNA_processed <- paste0(dir_resources, "snRNA_Processed_Data/")
dir_scRNA_auto <- paste0(dir_snRNA_processed, "scRNA_auto/")
dir_scRNA_auto_out <- paste0(dir_scRNA_auto, "outputs/")

# create output directory ----------------------------------------------------------
dir_analysis_results1 <- paste0(dir_analysis_results, "integration/")
dir.create(dir_analysis_results1)
dir_analysis_results2 <- paste0(dir_analysis_results1, "30_aliquot_integration/")
dir.create(dir_analysis_results2)
dir_analysis_results3 <- paste0(dir_analysis_results2, "docker_run_integration/")
dir.create(dir_analysis_results3)
dir_out <- paste0(dir_analysis_results3, run_id, "/")
dir.create(dir_out)

# input seurat processing summary to get the objects to process------------------------------------------------
path_seurat_summary <- paste0(dir_snRNA_processed, "scRNA_auto/summary/")
seurat_summary <- fread(input = paste0(path_seurat_summary,"ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200207.v1.tsv"), data.table = F)
## filter for the to be processed aliquots
seurat_summary2process <- seurat_summary %>%
  filter(Proceed_for_downstream == "Yes") %>%
  mutate(Path_seurat_object = paste0(dir_scRNA_auto_out, Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object

# input the seurat object and store in a list--------------------------------------------------------
renal.list <- list()
for (i in 1:nrow(seurat_summary2process)) {
  sample_id_tmp <- seurat_summary2process$Aliquot[i]
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[i]
  seurat_obj <- readRDS(file = seurat_obj_path)
  renal.list[[sample_id_tmp]] <- seurat_obj
}
length(renal.list)
rm(seurat_obj)

# Run the standard workflow for visualization and clustering ------------
### normalization
renal.list <- lapply(X = renal.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
## integrate without anchor
renal.integrated <- merge(x = renal.list[[1]], y = renal.list[2:length(renal.list)], project = "integrated")
rm(renal.list)
## get variably expressed genes
renal.integrated <- FindVariableFeatures(renal.integrated, selection.method = "vst", nfeatures = 3000)
## scale data with all the features
renal.integrated <- ScaleData(renal.integrated, features = rownames(renal.integrated@assays$RNA@counts)) 
## keep it consistant with individual processing pipeline
renal.integrated <- RunPCA(renal.integrated, npcs = 50, verbose = FALSE)
renal.integrated <- RunUMAP(renal.integrated, reduction = "pca", dims = 1:30)
renal.integrated <- FindNeighbors(renal.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
renal.integrated <- FindClusters(renal.integrated, resolution = 0.5)
## save as RDS file
saveRDS(object = renal.integrated, file = paste0(dir_out, "30_aliquot_integration.", run_id, ".RDS"), compress = T)

# plot dimensions by cluster with all aliquots together ------------------------------
## make sure the grouping variable is in the meta data

p <- DimPlot(renal.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)

file2write <- paste(dir_out, "Dimplot_by_Clusters.", run_id, ".png", sep="")
png(file = file2write, width = 1000, height = 800, res = 150)
print(p)
dev.off()

# plot dimensions by aliquot with all aliquots together ------------------------------
p <- DimPlot(renal.integrated, reduction = "umap", group.by = "orig.ident", label = F)

file2write <- paste(dir_out, "Dimplot_by_Aliquots.", run_id, ".png", sep="")
png(file = file2write, width = 1100, height = 800, res = 150)
print(p)
dev.off()

# plot dimensions by aliquot and by cluster ------------------------------
p <- DimPlot(renal.integrated, reduction = "umap", label = T, split.by = "orig.ident", ncol = 8)

file2write <- paste(dir_out, "Dimplot_by_Clusters_and_Aliquots.", run_id, ".png", sep="")
png(file = file2write, width = 5000, height = 2000, res = 150)
print(p)
dev.off()
