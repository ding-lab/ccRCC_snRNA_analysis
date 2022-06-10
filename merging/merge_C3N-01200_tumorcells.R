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
  "Seurat"
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

# preprocess --------------------------------------------------------------
paths_rds_df <- data.frame(sample_id = paste0("C3N-01200-", c("T1", "T2", "T3")),
                           path = paste0("./Data_Freezes/V2/snRNA/Tumor_Cell_Reclustered/C3N-01200-", c("T1", "T2", "T3"), ".tumorcellreclustered.20201124.v1.RDS"))

# process -----------------------------------------------------------------
## scale data with all the features
renal.list <- list()
for (i in 1:nrow(paths_rds_df)) {
  sample_id_tmp <- paths_rds_df$sample_id[i]
  seurat_obj_path <- paths_rds_df$path[i]
  seurat_obj <- readRDS(file = seurat_obj_path)
  ## take out the doublets
  barcode2scrublet_tmp_df <- barcode2scrublet_df %>%
    filter(Aliquot == sample_id_tmp) %>%
    filter(Barcode %in% rownames(seurat_obj@meta.data)) %>%
    filter(predicted_doublet)
  barcodes_scrublet <- barcode2scrublet_tmp_df$Barcode
  barcodes_keep <- rownames(seurat_obj@meta.data)
  barcodes_keep <- barcodes_keep[!(barcodes_keep %in% barcodes_scrublet)]
  ###
  print("subsetting")
  seurat_sub_obj <- subset(x = seurat_obj, cells = barcodes_keep)
  print(dim(seurat_sub_obj))
  renal.list[[sample_id_tmp]] <- seurat_sub_obj
}
length(renal.list)
rm(seurat_obj)

renal.list <- lapply(X = renal.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
## integrate without anchor
renal.integrated <- merge(x = renal.list[[1]], y = renal.list[2:length(renal.list)], project = "merged")
rm(renal.list)

renal.integrated <- SCTransform(renal.integrated, vars.to.regress = c("nCount_RNA","percent.mito"), return.only.var.genes = F)
cat("Finished SCTransform!\n")

## keep it consistant with individual processing pipeline
# renal.integrated <- RunPCA(renal.integrated, npcs = 30, verbose = FALSE)
renal.integrated <- RunPCA(renal.integrated, npcs = 45, verbose = FALSE)
cat("Finished RUNPCA!\n")

# determine the number of PCs --------------------------------------------------------------
p <- ElbowPlot(renal.integrated, ndims = 45)
file2write <- paste0(dir_out, "C3N-01200.Tumorcells.Merged.", run_id, ".pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
print(p)
dev.off()

file2write <- paste0(dir_out, "C3N-01200.Tumorcells.Merged.", run_id, ".txt")
sink(file2write)
# Determine percent of variation associated with each PC
pct <- renal.integrated[["pca"]]@stdev / sum(renal.integrated[["pca"]]@stdev) * 100
print(pct)

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
print(cumu)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
print(co1)
sink()


# finish processing -------------------------------------------------------
renal.integrated <- RunUMAP(renal.integrated, reduction = "pca", dims = 1:39)
cat("Finished RUNUMAP!\n")
renal.integrated <- FindNeighbors(renal.integrated, reduction = "pca", dims = 1:39, force.recalc = T)
cat("Finished FindNeighbors!\n")
renal.integrated <- FindClusters(renal.integrated, resolution = 0.5)
cat("Finished FindClusters!\n")
## save as RDS file
path_final_file <- paste0(dir_out, "C3N-01200.Tumorcells.Merged.", run_id, ".RDS")
saveRDS(object = renal.integrated, file = path_final_file, compress = T)
cat("Finished saving the output!\n")

# plot UMAP and fetch data ------------------------------------------------
Idents(renal.integrated) <- "orig.ident"
p <- DimPlot(renal.integrated)
file2write <- paste0(dir_out, "C3N-01200.Tumorcells.Merged.UMAP_bysample.", run_id, ".pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
print(p)
dev.off()
metadata_df <- FetchData(object = renal.integrated, vars = c("UMAP_1", "UMAP_2", "orig.ident"))
metadata_df$barcode_merged <- rownames(metadata_df)
file2write <- paste0(dir_out, "C3N-01200.Tumorcells.Merged.UMAP_data.", run_id, ".tsv")
write.table(x = metadata_df, file = file2write, quote = F, sep = "\t", row.names = F)
