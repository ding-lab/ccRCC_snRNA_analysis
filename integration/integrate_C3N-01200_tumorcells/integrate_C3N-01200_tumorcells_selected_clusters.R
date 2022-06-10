# Yige Wu @WashU Mar 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input doublet prediction
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)

# preprocess --------------------------------------------------------------
paths_rds_df <- data.frame(sample_id = paste0("C3N-01200-", c("T1", "T2", "T3")),
                           path = paste0("./Data_Freezes/V2/snRNA/Tumor_Cell_Reclustered/C3N-01200-", c("T1", "T2", "T3"), ".tumorcellreclustered.20201124.v1.RDS"))
unique(barcode2tumorsubcluster_df$Cluster_Name)
clusters_selected <- paste0("C3N-01200-", c("T1_C1", "T1_C2", "T1_C3", "T1_C4",
                                            "T2_C1", "T2_C2",
                                            "T3_C1"))

# process -----------------------------------------------------------------
## scale data with all the features
srat_list <- list()
for (i in 1:nrow(paths_rds_df)) {
  sample_id_tmp <- paths_rds_df$sample_id[i]
  seurat_obj_path <- paths_rds_df$path[i]
  seurat_obj <- readRDS(file = seurat_obj_path)
  ## take out the doublets
  barcode2scrublet_tmp_df <- barcode2scrublet_df %>%
    filter(Aliquot_WU == sample_id_tmp) %>%
    filter(predicted_doublet == T)
  barcodes_scrublet <- barcode2scrublet_tmp_df$Barcode
  barcodes_keep <- rownames(seurat_obj@meta.data); print(length(barcodes_keep))
  barcodes_keep <- barcodes_keep[!(barcodes_keep %in% barcodes_scrublet)]; print(length(barcodes_keep))
  barcodes_keep <- barcodes_keep[barcodes_keep %in% barcode2tumorsubcluster_df$barcode[barcode2tumorsubcluster_df$easy_id == sample_id_tmp & barcode2tumorsubcluster_df$Cluster_Name %in% clusters_selected]]; print(length(barcodes_keep))
  ###
  print("subsetting")
  seurat_sub_obj <- subset(x = seurat_obj, cells = barcodes_keep)
  print(dim(seurat_sub_obj))
  srat_list[[sample_id_tmp]] <- seurat_sub_obj
}
length(srat_list)
rm(seurat_obj)

# Run the standard workflow for visualization and clustering ------------
srat_list <- lapply(X = srat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
cat("Finished NormalizeData and FindVariableFeatures for the list!\n")

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = srat_list)
cat("Finished SelectIntegrationFeatures!\n")

srat_list <- lapply(X = srat_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

## chose one male and one female as references to scale up the integration
srat_anchors <- FindIntegrationAnchors(object.list = srat_list, anchor.features = features, reduction = "rpca", dims = 1:40, reference = 1)
cat("Finished FindIntegrationAnchors!\n")

# this command creates an 'integrated' data assay
srat_integrated <- IntegrateData(anchorset = srat_anchors, dims = 1:40)
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

# determine the number of PCs --------------------------------------------------------------
p <- ElbowPlot(srat_integrated, ndims = 45)
file2write <- paste0(dir_out, "C3N-01200.Tumorcells.Integrated.", run_id, ".pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
print(p)
dev.off()

file2write <- paste0(dir_out, "C3N-01200.Tumorcells.Integrated.", run_id, ".txt")
sink(file2write)
# Determine percent of variation associated with each PC
pct <- srat_integrated[["pca"]]@stdev / sum(srat_integrated[["pca"]]@stdev) * 100
print(pct)

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
print(cumu)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
print(co1)
sink()


# finish processing -------------------------------------------------------
srat_integrated <- RunUMAP(srat_integrated, reduction = "pca", dims = 1:39)
cat("Finished RUNUMAP!\n")
srat_integrated <- FindNeighbors(srat_integrated, reduction = "pca", dims = 1:39, force.recalc = T)
cat("Finished FindNeighbors!\n")
srat_integrated <- FindClusters(srat_integrated, resolution = 0.5)
cat("Finished FindClusters!\n")
## save as RDS file
path_final_file <- paste0(dir_out, "C3N-01200.Tumorcells.Integrated.", run_id, ".RDS")
saveRDS(object = srat_integrated, file = path_final_file, compress = T)
cat("Finished saving the output!\n")

# plot UMAP and fetch data ------------------------------------------------
Idents(srat_integrated) <- "orig.ident"
p <- DimPlot(srat_integrated)
file2write <- paste0(dir_out, "C3N-01200.Tumorcells.Integrated.UMAP_bysample.", run_id, ".pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
print(p)
dev.off()
metadata_df <- FetchData(object = srat_integrated, vars = c("UMAP_1", "UMAP_2", "orig.ident"))
metadata_df$barcode_merged <- rownames(metadata_df)
file2write <- paste0(dir_out, "C3N-01200.Tumorcells.Integrated.UMAP_data.", run_id, ".tsv")
write.table(x = metadata_df, file = file2write, quote = F, sep = "\t", row.names = F)
