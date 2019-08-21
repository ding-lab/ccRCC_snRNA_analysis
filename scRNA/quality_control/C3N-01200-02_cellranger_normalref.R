# Yige Wu @ WashU 20189 Aug
## for QC the snRNA-seq data for C3N-01200-02


# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# Set IDs -----------------------------------------------------------------
sample_id <- "C3N-01200-02_normalref"

# Set path ----------------------------------------------------------------
matrix_dir <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/raw_feature_bc_matrix/")
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
## get ensemble gene ids to gene symbol mapping
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
rownames(feature.names) <- feature.names$V1
colnames(feature.names) <- c("ensemble_gene_id", "symbol", "assay")
feature.names %>% head()
feature.names %>% nrow()

# Read input --------------------------------------------------------------
## readMM: Read matrices stored in the Harwell-Boeing or MatrixMarket formats or write sparseMatrix objects to one of these formats.
input <- readMM(file = matrix.path)

## get cell barcodes
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names %>% head
barcode.names %>% nrow()

## make the column names of the expression matrix to be each cell barcode
## make the row names of the expression matrix to be each gene
colnames(input) = barcode.names$V1
rownames(input) = feature.names$ensemble_gene_id

###########################################
######## ANALYSIS
###########################################
# calculate the expression count sum of each cell -------------------------
bc_sums <- colSums(input)
summary(bc_sums)

## filter by expression sum and get the cell barcodes with expression sum over 300 (#cells: 77870)
bg_id <- names(bc_sums[bc_sums >= 300])
bg_id %>% length()

# input = input[,bg_id]
# saveRDS(object = input, file = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_matrix_counts_over300_per_cell.RDS"))

## filter by expression sum and get the cell barcodes with expression sum over 1000 under 50000 (#cells: 12715)
bg_id <- names(bc_sums[bc_sums >= 1000 & bc_sums <= 50000])
bg_id %>% length()

input = input[,bg_id]
saveRDS(object = input, file = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_matrix_counts_over1000_under50000_per_cell.RDS"))

## filter by expression sum and get the cell barcodes with expression sum over 1500 under 50000 (#cells: 7127)
bg_id <- names(bc_sums[bc_sums >= 1500 & bc_sums <= 50000])
bg_id %>% length()

input = input[,bg_id]
saveRDS(object = input, file = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_matrix_ccounts_over1500_under50000_per_cell.RDS"))

## filter by expression sum and get the cell barcodes with expression sum over 2000 under 50000 (#cells: 4727)
bg_id <- names(bc_sums[bc_sums >= 2000 & bc_sums <= 50000])
bg_id %>% length()

# input = input[,bg_id]
# saveRDS(object = input, file = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_matrix_ccounts_over2000_under50000_per_cell.RDS"))

rm(bc_sums)
rm(bg_id)

# if run after the counts are filtered --------------------------------------------------
# file2input <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_matrix_ccounts_over1500_under50000_per_cell.RDS")
file2input <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_matrix_counts_over1000_under50000_per_cell.RDS")
input <- readRDS(file = file2input)
dim(input)
primary_object <- CreateSeuratObject(counts = input,project="Kidney", min.cells = 0)

# calculate the the percentage of mitochondrial gene expression ---------------
primary_object.row.symbol <- feature.names[match(rownames(primary_object@assays$RNA@data), rownames(input)), 2]
mito.genes <- grep(pattern = "^MT-", x = primary_object.row.symbol, value=F)
mito.genes
percent.mito <- Matrix::colSums(primary_object@assays$RNA@data[mito.genes, ])/Matrix::colSums(primary_object@assays$RNA@data)
primary_object <- AddMetaData(object = primary_object, metadata = percent.mito, col.name = "percent.mito")
# file2write <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_primary_seurat_object_counts_over1500_under50000_per_cell.RDS")
file2write <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_primary_seurat_object_counts_over1000_under50000_per_cell.RDS")
saveRDS(object = primary_object, file = file2write)

# plot Pre-filtering QC metrics ---------------------------------------------------------
# file2write <- paste(makeOutDir(), "PreQC_in_sample_", sample_id, ".pdf", sep="")
# pdf(file = file2write, width=15, height=9)
# VlnPlot(object = primary_object, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3)
# dev.off()

# # get cutoffs
# info <- as.data.frame(x=FetchData(primary_object, vars.all=c("nCount_RNA","nFeature_RNA", "percent.mito")))
# 
# m <- vector()
# n <- vector()
# for(i in seq(300, 3000, 100))
# {
#   tmp <- info[info$nUMI > i, ]
#   m <- c(m, sum(tmp$percent.mito > 0.1)/nrow(tmp))
#   n <- c(n, nrow(tmp))
# }
# data <- data.frame(cutoff=seq(300,3000, 100), percentage=m, number=n)
# 
# pdf(paste("Different_UMI_cutoff_in_sample_", sample_id, ".pdf", sep=""))
# par(mar = c(5, 5, 3, 5))
# plot(data$cutoff, data$percentage, type ="l", ylab = "Percentage of dead cells",
#      main = "", xlab = "UMI cutoff",
#      col = "blue", las=1)
# par(new = TRUE)
# plot(data$cutoff, data$number, type = "l", xaxt = "n", yaxt = "n",
#      ylab = "", xlab = "", col = "red", lty = 2)
# axis(side = 4)
# mtext("Number of cells", side = 4, line = 3)
# legend("topright", c("Percentage", "Number"), col = c("blue", "red"), lty = c(1, 2))
# dev.off()


# filtering primary seurat object ------------------------------------------------
qced_object <- subset(primary_object, subset = nCount_RNA > 1000 &  nCount_RNA < 50000)
qced_object@meta.data %>%
  dim()

qced_object <- subset(primary_object, subset = nCount_RNA > 1500 &  nCount_RNA < 50000)
qced_object@meta.data %>%
  dim()

qced_object <- subset(primary_object, subset = nCount_RNA > 1000 &  nCount_RNA < 50000 & percent.mito < 0.10)
qced_object@meta.data %>%
  dim()
# [1] 3708    4
qced_object <- subset(primary_object, subset = nCount_RNA > 1500 &  nCount_RNA < 50000 & percent.mito < 0.10)
qced_object@meta.data %>%
  dim()

qced_object <- subset(primary_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 1000 &  nCount_RNA < 50000 & percent.mito < 0.10)
qced_object@meta.data %>%
  dim()
# [1] 3703    4
qced_object <- subset(primary_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 1500 &  nCount_RNA < 50000 & percent.mito < 0.10)
qced_object@meta.data %>%
  dim()
# [1] 1475    4
# file2write <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_QCed_seurat_object_nFeature_over200_under5000_nCount_over1500_under50000_percent.mito_0.1.RDS")
# saveRDS(object = qced_object, file = file2write)
file2write <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_QCed_seurat_object_nFeature_over200_under5000_nCount_over1000_under50000_percent.mito_0.1.RDS")
saveRDS(object = qced_object, file = file2write)

# plot after filtering QC metrics -----------------------------------------
file2write <- paste(makeOutDir(), "Post-Filtering_QC_metrics_in_sample_", sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=15, height=9)
VlnPlot(object = qced_object, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3)
dev.off()

# ##Section 3: Downstream analysis ----------------------------------------
# normalization
qced_object <- NormalizeData(qced_object, normalization.method = "LogNormalize", scale.factor = 10000)

# variable genes
qced_object <- FindVariableFeatures(object = qced_object, selection.method = "vst", nfeatures = 2000)

# scale
qced_object <- SCTransform(qced_object, vars.to.regress = c("nCount_RNA","percent.mito")) 

# PCA
qced_object <- RunPCA(qced_object)
# file2write <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_PCAed_seurat_object_nFeature_over200_under5000_nCount_over1500_under50000_percent.mito_0.1.RDS")
# saveRDS(object = qced_object, file = file2write)
# file2read <- file2write
# qced_object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/PCAed_Seurat_object.RDS")

## run umap
## pip install umap-learn
qced_object <- RunUMAP(qced_object, dims= 1:30)
# saveRDS(object = qced_object, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/UMAPed_Seurat_object.RDS")

## find neighbors
qced_object <- FindNeighbors(qced_object, dims = 1:30)

## find clusters
qced_object <- FindClusters(qced_object)
# file2write <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_Clustered_seurat_object_nFeature_over200_under5000_nCount_over1500_under50000_percent.mito_0.1.RDS")
file2write <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_Clustered_seurat_object_nFeature_over200_under5000_nCount_over1000_under50000_percent.mito_0.1.RDS")
saveRDS(object = qced_object, file = file2write)

# write out the expression matrics for pollock
fwrite(x = as.data.frame(qced_object@assays$SCT@counts), row.names = TRUE, file = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/", sample_id, "/", sample_id, "_formatted_raw.tsv"),sep="\t")

# Plot dims ---------------------------------------------------------------
file2write <- paste(makeOutDir(), "DimPlot_raw_cell_clusters_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=12, height=12)
DimPlot(qced_object, label = TRUE) + NoLegend() + theme_bw()
dev.off()

# Tumor marker plot -------------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_Tumor_marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=12, height=15)
top.markers <- c("NPHS2", "LRP2", "AQP1", "AQP2", "UMOD")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Endomthelial marker plot -------------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_Endomthelial_marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=12, height=6)
top.markers <- c("VWF", "PECAM1", "FLT4", "FLT1", "FLT3", "KDR")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Tumor oncogene plot -----------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_Tumor_Oncogenes_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=15, height=12)
top.markers <- c("HIF1A", "EPAS1" , "MTOR", "MET", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
tmp@meta.data[,top.markers] <- GetAssayData(tmp)[c(markers_ens),] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Tumor TSG plot -----------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_Tumor_Suppressor_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=12, height=15)
top.markers <- c("VHL", "PBRM1" , "BAP1", "SETD2", "PTEN", "TP53")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
tmp@meta.data[,top.markers] <- GetAssayData(tmp)[c(markers_ens),] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Proximal Tubule marker plot -----------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_Proximal_Tubule_Marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=12, height=15)
top.markers <- c("SLC22A8", "SLC17A3" , "SLC22A7", "SLC16A9", "SLC7A13", "SLC34A1", "SLC13A3")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# CD4 marker plot ---------------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_CD4Tcell_Marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=15, height=12)
top.markers <- as.vector(CD4T_markers$gene_symbol)
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# CD8 marker plot ---------------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_CD8Tcell_Marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=12, height=15)
top.markers <- as.vector(CD8T_markers$gene_symbol)
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Macrophage marker plot ---------------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_Macrophage_Marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=12, height=12)
top.markers <- as.vector(Macrophage_markers$gene_symbol)
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Monocyte marker plot ---------------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_Monocyte_Marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=12, height=6)
top.markers <- as.vector(Monocyte_markers$gene_symbol)
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# NK marker plot ---------------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_NK_Marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=12, height=12)
top.markers <- as.vector(NK_markers$gene_symbol)
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Plasma marker plot ---------------------------------------------------------
file2write <- paste(makeOutDir(), "FeaturePlot_Plasma_Marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
pdf(file = file2write, width=12, height=12)
top.markers <- as.vector(Plasma_markers$gene_symbol)
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Erythrocytes marker plot ---------------------------------------------------------
# file2write <- paste(makeOutDir(), "FeaturePlot_Erythrocytes_Marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
# pdf(file = file2write, width=12, height=12)
# top.markers <- as.vector(Erythrocytes_markers$gene_symbol)
# markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
# tmp <- qced_object
# # Save assay data to metadata
# markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
# tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
# FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
# dev.off()


# DC marker plot ---------------------------------------------------------
# file2write <- paste(makeOutDir(), "FeaturePlot_DC_Marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
# pdf(file = file2write, width=12, height=12)
# top.markers <- as.vector(DC_markers$gene_symbol)
# markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
# tmp <- qced_object
# # Save assay data to metadata
# markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
# tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
# FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
# dev.off()

# B marker plot ---------------------------------------------------------
# file2write <- paste(makeOutDir(), "FeaturePlot_B_Marker_expression_in_sample_",sample_id, ".20190816_v2.pdf", sep="")
# pdf(file = file2write, width=12, height=12)
# top.markers <- as.vector(B_markers$gene_symbol)
# markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
# tmp <- qced_object
# # Save assay data to metadata
# markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
# tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
# FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
# dev.off()


# # Section 4: Label cells --------------------------------------------------
# current_labels <- Idents(renal) %>% levels
# #!! REPLACE MANNUAL ASSIGNED CELL LABEL BELOW !!
# new_labels     <- c('CD8Tcell','CD4Tcell','Macrophage','NK(1)','Epithelial(rc, stem)','NK(2)','DC','NK(3)','Tumor(Megalin)','Treg','NK(4)','Endothelial','Astrocyte','Tumor(Resistance)','Bcell','Pro-Bcell','Unknown')
# Idents(renal)  <- plyr::mapvalues(as.character(Idents(renal)), from= current_labels ,to = new_labels)
# renal          <- AddMetaData(renal, Idents(renal), col.name = 'cell_type')  
# 
# # Show summary
# table(Idents(renal)) %>% as.data.frame()
# 
# 
# plot pollock output -----------------------------------------------------
celltypes <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02_normalref/C3N-01200-02_normalref_cell_types.tsv", data.table = F)
qced_object$predicted_cell_type = merge(as.data.frame(rownames(qced_object@meta.data)), celltypes, by.x = 'rownames(qced_object@meta.data)', by.y='id')$predicted_cell_type
file2write <- paste(makeOutDir(), "DimPlot_pollock_cell_type_in_sample_",sample_id, ".20190816_v1.pdf", sep="")
pdf(file = file2write, width=12, height=12)
DimPlot(qced_object, reduction = "umap", group.by = "predicted_cell_type", label=TRUE)
dev.off()

