# Yige Wu @ WashU 20189 Aug
## for QC the snRNA-seq data for C3N-01200-02


# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# Set IDs -----------------------------------------------------------------
sample_id <- "C3N-01200-02"

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

###########################################
######## FUNCTIONS
###########################################

filter_func <- function(gbm, id_to_remove)
{
        new_matrix <- exprs(gbm)[,!colnames(exprs(gbm)) %in% id_to_remove]
        pd <- data.frame(barcode = colnames(new_matrix))
        rownames(pd) <- colnames(new_matrix)
        return(newGeneBCMatrix(new_matrix, fd = fData(gbm), pd = pd, template = gbm))
}



run_scale_pca_cluster = function( obj, name ) {
  # takes a seurat object with raw data and perform normalization, PCA, clustering and tSNE
  # obj = seurat object
  # name = file name prefix
  # returns modified seurat object
  
  # normalization
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # variable genes
  pdf(paste("Variable_genes_in_sample_", name, "_", sample_id, ".pdf", sep=""))
  obj <- FindVariableGenes(object = obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = T)
  dev.off()
  length(x = obj@var.genes)
  
  # scale
  obj <- ScaleData(obj, vars.to.regress = c("nUMI", "percent.mito"))
  
  # PCA
  obj <- RunPCA(obj, pc.genes = obj@var.genes, do.print = F)
  
  # cluster
  obj <- FindClusters(obj, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
  
  # tSNE
  obj <- RunTSNE(obj, dims.use = 1:10, do.fast = TRUE)
  pdf(paste("TSNE_in_sample_", name, "_", sample_id, ".pdf", sep=""))
  TSNEPlot(obj, pt.size=2, label.size = 6)
  dev.off()
  
  # DEG (marker)
  obj.markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  obj.markers$gene <- feature.names[match(obj.markers$gene, feature.names$id), 2]
  obj.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
  DEdata = obj.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
  write.table(DEdata, file = "./", name, "_", sample_id, "DEdata.txt", quote = F, sep = "\t")
  
  # Heatmap DEG
  pdf(paste("Heatmap_DEG_in_sample_", name, "_", sample_id, ".pdf", sep=""))
  top.markers <- primary_object.markers %>% group_by(cluster) %>% top_n(4, avg_logFC)
  top.markers = top.markers$gene
  markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
  tmp <- primary_object
  tmp@scale.data <- tmp@scale.data[rownames(tmp@scale.data) %in% markers_ens, ]
  rownames(tmp@scale.data) <- feature.names[match(rownames(tmp@scale.data), feature.names$id), 2]
  # setting slim.col.label to TRUE will print just the cluster IDS instead of
  # every cell name
  DoHeatmap(object = tmp, genes.use = top.markers, slim.col.label = TRUE, remove.key = TRUE, group.cex = 8)
  dev.off()
  
  return( obj )
}

###########################################
######## ANALYSIS
###########################################
# library -----------------------------------------------------------------
# source("./Box/Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")
# 
# 
# # Read input --------------------------------------------------------------
# ## readMM: Read matrices stored in the Harwell-Boeing or MatrixMarket formats or write sparseMatrix objects to one of these formats.
# input <- readMM(file = matrix.path)
# 
# ## get cell barcodes
# barcode.names = read.delim(barcode.path,
#                            header = FALSE,
#                            stringsAsFactors = FALSE)
# barcode.names %>% head
# barcode.names %>% nrow()
# 
# ## make the column names of the expression matrix to be each cell barcode
# ## make the row names of the expression matrix to be each gene
# colnames(input) = barcode.names$V1
# rownames(input) = feature.names$V1
# 
# ## calculate the expression count sum of each cell
# bc_sums <- colSums(input)
# summary(bc_sums)
# 
# ## filter by expression sum and get the cell barcodes with expression sum over 300 (#cells: 81444)
# bg_id <- names(bc_sums[bc_sums >= 300])
# bg_id %>% length()
# 
# input = input[,bg_id]
# saveRDS(object = input, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/counts_over300_per_cell.RDS")
# 
# ## filter by expression sum and get the cell barcodes with expression sum over 1000 under 50000 (#cells: 78575)
# bg_id <- names(bc_sums[bc_sums >= 1000 & bc_sums <= 50000])
# bg_id %>% length()
# 
# input = input[,bg_id]
# saveRDS(object = input, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/counts_over1000_under50000_per_cell.RDS")
# 
# ## filter by expression sum and get the cell barcodes with expression sum over 1500 under 50000 (#cells: 43987)
# bg_id <- names(bc_sums[bc_sums >= 1500 & bc_sums <= 50000])
# bg_id %>% length()
# 
# input = input[,bg_id]
# saveRDS(object = input, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/counts_over1500_under50000_per_cell.RDS")
# 
# ## filter by expression sum and get the cell barcodes with expression sum over 2000 under 50000 (#cells: 23056)
# bg_id <- names(bc_sums[bc_sums >= 2000 & bc_sums <= 50000])
# bg_id %>% length()
# 
# input = input[,bg_id]
# saveRDS(object = input, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/counts_over2000_under50000_per_cell.RDS")
# 
# rm(bc_sums)
# rm(bg_id)

# if run after the counts are filtered --------------------------------------------------
# input <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/counts_over300_per_cell.RDS")
## Error: vector memory exhausted (limit reached?)
# input <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/counts_over1000_under50000_per_cell.RDS")
## Error: vector memory exhausted (limit reached?)
input <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/counts_over1500_under50000_per_cell.RDS")
dim(input)
primary_object <- CreateSeuratObject(counts = input,project="Kidney", min.cells = 0)

# calculate the the percentage of mitochondrial gene expression ---------------
primary_object.row.symbol <- feature.names[match(rownames(primary_object@assays$RNA@data), rownames(input)), 2]
mito.genes <- grep(pattern = "^MT-", x = primary_object.row.symbol, value=F)
percent.mito <- Matrix::colSums(primary_object@assays$RNA@data[mito.genes, ])/Matrix::colSums(primary_object@assays$RNA@data)
primary_object <- AddMetaData(object = primary_object, metadata = percent.mito, col.name = "percent.mito")
# saveRDS(object = primary_object, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/Primary_Seurat_object.RDS")

# plot QC metrics ---------------------------------------------------------
# fp_tmp <- paste(makeOutDir(), "QC_in_sample_",sample_id, ".pdf", sep="")
# pdf(file = fp_tmp, width=15, height=9)
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


# QC primary seurat object ------------------------------------------------
# qced_object <- subset(primary_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 10000 &  nCount_RNA < 50000 & percent.mito < 0.05)
# qced_object@meta.data %>%
#   dim()
# qced_object <- subset(primary_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 5000 &  nCount_RNA < 50000 & percent.mito < 0.05)
# qced_object@meta.data %>%
#   dim()
qced_object <- subset(primary_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 4000 &  nCount_RNA < 50000 & percent.mito < 0.05)
qced_object@meta.data %>%
  dim()
# qced_object <- subset(primary_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 2500 &  nCount_RNA < 50000 & percent.mito < 0.05)
# qced_object@meta.data %>%
#   dim()
# saveRDS(object = qced_object, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/QCed_Seurat_object.RDS")
rm(primary_object)


# plot after filtering QC metrics -----------------------------------------
# fp_tmp <- paste(makeOutDir(), "After_QC_in_sample_",sample_id, ".20190815_v2.pdf", sep="")
# pdf(file = fp_tmp, width=15, height=9)
# VlnPlot(object = qced_object, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3)
# dev.off()


# ##Section 3: Downstream analysis ----------------------------------------
# normalization
qced_object <- NormalizeData(qced_object, normalization.method = "LogNormalize", scale.factor = 10000)

# variable genes
qced_object <- FindVariableFeatures(object = qced_object, selection.method = "vst", nfeatures = 2000)

# scale
qced_object <- SCTransform(qced_object, vars.to.regress = c("nCount_RNA","percent.mito")) 

# PCA
qced_object <- RunPCA(qced_object)
# saveRDS(object = qced_object, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/PCAed_Seurat_object.RDS")
qced_object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/PCAed_Seurat_object.RDS")

## run umap
## pip install umap-learn
qced_object <- RunUMAP(qced_object, dims= 1:30)
# saveRDS(object = qced_object, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/UMAPed_Seurat_object.RDS")

## find neighbors
qced_object <- FindNeighbors(qced_object, dims = 1:30)

## find clusters
qced_object <- FindClusters(qced_object)
# saveRDS(object = qced_object, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/Clustered_Seurat_object.RDS")

## write out the expression matrics for pollock
# fwrite(x = as.data.frame(qced_object@assays$SCT@counts), row.names = TRUE, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02/C3N-01200-02_formatted_raw.tsv",sep="\t")

# Plot dims ---------------------------------------------------------------
fp_tmp <- paste(makeOutDir(), "DimPlot_raw_cell_clusters_in_sample_",sample_id, ".20190815_v1.pdf", sep="")
pdf(file = fp_tmp, width=12, height=12)
DimPlot(qced_object, label = TRUE) + NoLegend() + theme_bw()
dev.off()

fp_tmp <- paste(makeOutDir(), "DimPlot_CD14_CD16_marker_expression_in_sample_",sample_id, "_ens_id.20190815_v2.pdf", sep="")
pdf(file = fp_tmp, width=12, height=6)
top.markers <- c("CD14", "FCGR3A", "LYZ")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Tumor marker plot -------------------------------------------------------
fp_tmp <- paste(makeOutDir(), "DimPlot_Tumor_marker_expression_in_sample_",sample_id, "_ens_id.20190815_v1.pdf", sep="")
pdf(file = fp_tmp, width=12, height=12)
top.markers <- c("NPHS2", "LRP2", "AQP1", "AQP2", "UMOD")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
FeaturePlot(object = tmp, features = markers_ens, cols = c("grey", "red"), reduction = "umap")
dev.off()

fp_tmp <- paste(makeOutDir(), "DimPlot_Tumor_marker_expression_in_sample_",sample_id, "_ens_id.20190815_v2.pdf", sep="")
pdf(file = fp_tmp, width=12, height=6)
top.markers <- c("NPHS2", "LRP2", "AQP1", "AQP2", "UMOD")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Tumor oncogene plot -----------------------------------------------------
#Create temprary Seurat object and assign
fp_tmp <- paste(makeOutDir(), "DimPlot_Tumor_Oncogenes_expression_in_sample_",sample_id, "_ens_id.20190815_v1.pdf", sep="")
pdf(file = fp_tmp, width=12, height=12)
top.markers <- c("HIF1A", "EPAS1" , "MTOR", "MET", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
FeaturePlot(object = tmp, features = markers_ens, cols = c("grey", "red"), reduction = "umap")
dev.off()

fp_tmp <- paste(makeOutDir(), "DimPlot_Tumor_Oncogenes_expression_in_sample_",sample_id, "_ens_id.20190815_v2.pdf", sep="")
pdf(file = fp_tmp, width=15, height=12)
top.markers <- c("HIF1A", "EPAS1" , "MTOR", "MET", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
tmp@meta.data[,top.markers] <- GetAssayData(tmp)[c(markers_ens),] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Tumor TSG plot -----------------------------------------------------
fp_tmp <- paste(makeOutDir(), "DimPlot_Tumor_Suppressor_expression_in_sample_",sample_id, "_ens_id.20190815_v1.pdf", sep="")
pdf(file = fp_tmp, width=15, height=12)
top.markers <- c("VHL", "PBRM1" , "BAP1", "SETD2", "PTEN", "TP53")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
tmp@meta.data[,top.markers] <- GetAssayData(tmp)[c(markers_ens),] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# Proximal Tubule marker plot -----------------------------------------------------
fp_tmp <- paste(makeOutDir(), "DimPlot_Proximal_Tubule_Marker_expression_in_sample_",sample_id, "_ens_id.20190815_v1.pdf", sep="")
pdf(file = fp_tmp, width=15, height=12)
top.markers <- c("SLC22A8", "SLC17A3" , "SLC22A7", "SLC16A9", "SLC7A13", "SLC34A1", "SLC13A3")
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()


# CD4 marker plot ---------------------------------------------------------
fp_tmp <- paste(makeOutDir(), "DimPlot_CD4Tcell_Marker_expression_in_sample_",sample_id, ".20190815_v1.pdf", sep="")
pdf(file = fp_tmp, width=15, height=12)
top.markers <- as.vector(CD4T_markers$gene_symbol)
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# CD8 marker plot ---------------------------------------------------------
fp_tmp <- paste(makeOutDir(), "DimPlot_CD8Tcell_Marker_expression_in_sample_",sample_id, ".20190815_v1.pdf", sep="")
pdf(file = fp_tmp, width=15, height=12)
top.markers <- as.vector(CD8T_markers$gene_symbol)
markers_ens <- feature.names[match(top.markers, feature.names$symbol), 1]
tmp <- qced_object
# Save assay data to metadata
markers_ens2plot <- intersect(rownames(tmp@assays$SCT@meta.features), markers_ens)
tmp@meta.data[,feature.names[markers_ens2plot, "symbol"]] <- GetAssayData(tmp)[markers_ens2plot,] %>% t()
FeaturePlot(object = tmp, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()


# Section 4: Label cells --------------------------------------------------
current_labels <- Idents(renal) %>% levels
#!! REPLACE MANNUAL ASSIGNED CELL LABEL BELOW !!
new_labels     <- c('CD8Tcell','CD4Tcell','Macrophage','NK(1)','Epithelial(rc, stem)','NK(2)','DC','NK(3)','Tumor(Megalin)','Treg','NK(4)','Endothelial','Astrocyte','Tumor(Resistance)','Bcell','Pro-Bcell','Unknown')
Idents(renal)  <- plyr::mapvalues(as.character(Idents(renal)), from= current_labels ,to = new_labels)
renal          <- AddMetaData(renal, Idents(renal), col.name = 'cell_type')  

# Show summary
table(Idents(renal)) %>% as.data.frame()


# plot pollock output -----------------------------------------------------
celltypes <- as.data.frame(pollock_tab)
qced_object$predicted_cell_type = merge(as.data.frame(rownames(qced_object@meta.data)), celltypes, by.x = 'rownames(qced_object@meta.data)', by.y='id')$predicted_cell_type
fp_tmp <- paste(makeOutDir(), "DimPlot_pollock_cell_type_in_sample_",sample_id, ".20190815_v1.pdf", sep="")
pdf(file = fp_tmp, width=12, height=12)
DimPlot(qced_object, reduction = "umap", group.by = "predicted_cell_type", label=TRUE)
dev.off()
