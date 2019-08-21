library(cellrangerRkit)
library(Seurat)
library(dplyr)

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
  obj.markers$gene <- fData(input)[match(obj.markers$gene, fData(input)$id), 2]
  obj.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
  DEdata = obj.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
  write.table(DEdata, file = "./", name, "_", sample_id, "DEdata.txt", quote = F, sep = "\t")
  
  # Heatmap DEG
  pdf(paste("Heatmap_DEG_in_sample_", name, "_", sample_id, ".pdf", sep=""))
  top.markers <- panc.markers %>% group_by(cluster) %>% top_n(4, avg_logFC)
  top.markers = top.markers$gene
  markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
  tmp <- panc
  tmp@scale.data <- tmp@scale.data[rownames(tmp@scale.data) %in% markers_ens, ]
  rownames(tmp@scale.data) <- fData(input)[match(rownames(tmp@scale.data), fData(input)$id), 2]
  # setting slim.col.label to TRUE will print just the cluster IDS instead of
  # every cell name
  DoHeatmap(object = tmp, genes.use = top.markers, slim.col.label = TRUE, remove.key = TRUE, group.cex = 8)
  dev.off()
  
  return( obj )
}

###########################################
######## ANALYSIS
###########################################

# input
sample_id <- "TWAN-RCC-A1"
cell_path <- "/diskmnt/Datasets/Kidney/scRNA"
# cell_path = "~/DingLab/SingleCell/Data/Kidney_Data/"
input <- get_matrix_from_h5(paste(cell_path, "/", sample_id, "/outs/raw_gene_bc_matrices_h5.h5", sep=""))

# pre-filter
bc_sums <- colSums(input)
bg_id <- names(bc_sums[bc_sums < 300])


gbm <- filter_func(input, bg_id)
expr <- exprs(gbm)

# set up primary object
panc <- CreateSeuratObject(raw.data = expr, project = "Kidney", min.genes = 0, min.cells = 0)

# QC
panc.row.symbol <- fData(input)[match(rownames(panc@data), fData(input)$id), 2]
mito.genes <- grep(pattern = "^MT-", x = panc.row.symbol, value=F)
percent.mito <- Matrix::colSums(panc@raw.data[mito.genes, ])/Matrix::colSums(panc@raw.data)

panc <- AddMetaData(object = panc, metadata = percent.mito, col.name = "percent.mito")
pdf(paste("QC_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
VlnPlot(object = panc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# cutoff
info <- as.data.frame(x=FetchData(panc, vars.all=c("nGene","nUMI", "percent.mito")))

m <- vector()
n <- vector()
for(i in seq(300, 3000, 100))
{
  tmp <- info[info$nUMI > i, ]
  m <- c(m, sum(tmp$percent.mito > 0.1)/nrow(tmp))
  n <- c(n, nrow(tmp))
}
data <- data.frame(cutoff=seq(300,3000, 100), percentage=m, number=n)

pdf(paste("Different_UMI_cutoff_in_sample_", sample_id, ".pdf", sep=""))
par(mar = c(5, 5, 3, 5))
plot(data$cutoff, data$percentage, type ="l", ylab = "Percentage of dead cells",
     main = "", xlab = "UMI cutoff",
     col = "blue", las=1)
par(new = TRUE)
plot(data$cutoff, data$number, type = "l", xaxt = "n", yaxt = "n",
     ylab = "", xlab = "", col = "red", lty = 2)
axis(side = 4)
mtext("Number of cells", side = 4, line = 3)
legend("topright", c("Percentage", "Number"), col = c("blue", "red"), lty = c(1, 2))
dev.off()


# filter
panc <- FilterCells(object = panc, subset.names = c("nGene", "nUMI", "percent.mito"), low.thresholds = c(200, 1000, -Inf), high.thresholds = c(5000, 50000, 0.10))
nrow(FetchData(panc))
# TWAN-RCC-A1: 3623
pdf(paste("After_QC_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
VlnPlot(object = panc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# normalization
panc <- NormalizeData(panc, normalization.method = "LogNormalize", scale.factor = 10000)

# variable genes
pdf(paste("Variable_genes_in_sample_", sample_id, ".pdf", sep=""))
panc <- FindVariableGenes(object = panc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = T)
dev.off()
length(x = panc@var.genes)

# scale
panc <- ScaleData(panc, vars.to.regress = c("nUMI", "percent.mito"))

# PCA
panc <- RunPCA(panc, pc.genes = panc@var.genes, do.print = F)

# cluster
panc <- FindClusters(panc, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)

# tSNE
panc <- RunTSNE(panc, dims.use = 1:10, do.fast = TRUE)
pdf(paste("TSNE_in_sample_", sample_id, ".pdf", sep=""))
TSNEPlot(panc, pt.size=2, label.size = 6)
dev.off()

# DEG (marker)
panc.markers <- FindAllMarkers(object = panc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
panc.markers$gene <- fData(input)[match(panc.markers$gene, fData(input)$id), 2]
panc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
# cluster0 marker gene:IL7R(CD4+)
# cluster1 marker gene:GNLY, CD7(NK/lymphocytes)
# cluster2 marker gene:CD8A(CD8+ T cell)
# cluster3 marker gene:APOC1(Macrophages)
# cluster4 marker gene:CRYAB, MIF, KRT18(renal cell, high stem, epithelial)
# cluster5 marker gene:ITGAX(CD11c)(Macrophages+lymphocyte signature=DC)
# cluster6 marker gene:LRP2-megalin(Tumor)
# cluster7 marker gene:AQP1(Tumor+Endothelial-VWF)
# cluster8 marker gene:RGSS(Astrocytes)
# cluster9 marker gene:AQP1(Tumor+some resistance-EGR1, astrocyte-CLU)
# cluster10 marker gene:CD79A(B-cell)
# cluster11 marker gene:UBE2C(pro B-cell)
# cluster12 marker gene:COL1A1(Astrocytes/Mesangial cells)

#cluster0 marker gene:KRT19 (ductal)
#cluster1 marker gene:CD3D, CD8A (CD8+)
#cluster2 marker gene:HLA-DPA1 (MHC class II) 
#cluster3 marker gene:KRT19 (ductal) 
#cluster4 marker gene:CD3D, IL7R (CD4+)
#cluster5 marker gene:CD79A (B)
#cluster6 marker gene:KRT19 (ductal)
#cluster7 marker gene:TPSAB1 (mast)
#cluster8 marker gene:CD79A (B)
#cluster9 marker gene:COL1A1 (PSC)
#cluster10 marker gene:KRT19 (ductal)
#cluster11 marker gene:PRSS1 (acinar)
#cluster12 marker gene:PLVAP (Endothelial)

# CD3 - TCELL
# CD8A - 8+ T
# CD4, IL7R - CD4+
# CD79A, CD79B, BLK, MS4A1 (CD20) - BCELL
# CD14, CD163, CSF1R - Monocytes/Macrophage
# CD68, CD14, CD163 - Macrophage
# CD7 - NK
# FUT4 (CD15) - Granulocyte
# ITGAX (CD11c) - DC
# IL3RA (CD123) - pDC
# MRC1 (CD206) - Pro tumor macrophage
# NPHS2 (podocin), LRP2 (megalin), AQP1 & 2, UMOD (Uromodulin) - Tumor

# save
saveRDS(panc, file = paste("panc_backup_object_in_sample_", sample_id, ".rds", sep=""))

# marker plot
pdf(paste("Marker_expression_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("CD8A", "CD4", "IL7R", "CD3G", "CD79A", "CD79B", "BLK", "MS4A1", "CD14", "CD68", "CD163", "CSF1R", "CD7", "FUT4", "ITGAX", "IL3RA", "MRC1")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# tumor marker plot
pdf(paste("Tumor_marker_expression_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("NPHS2", "LRP2", "AQP1", "AQP2", "UMOD")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# stem plot
pdf(paste("Stem_marker_expression_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("PROM1", "CD44", "ENG")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()


# marker plot
pdf(paste("CellType_marker_expression_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("IL7R", "CD7", "CD8A", "APOC1", "CRYAB", "ITGAX", "LRP2", "AQP1", "CD79A", "UBE2C", "COL1A1")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# pdf(paste("Endocrine_marker_expression_in_sample_", sample_id, ".pdf", sep=""))
# second.markers <- c("GCG", "INS", "PPY", "SST", "GHRL")
# markers_ens <- fData(input)[match(second.markers, fData(input)$symbol), 1]
# tmp <- panc
# tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
# rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
# FeaturePlot(object = tmp, features.plot = second.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
# dev.off()
# 
# pdf(paste("Tumor_marker_expression_in_sample_", sample_id, ".pdf", sep=""))
# second.markers <- c("KRAS", "ERBB2", "GATA6", "MYC", "AKT2")
# markers_ens <- fData(input)[match(second.markers, fData(input)$symbol), 1]
# tmp <- panc
# tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
# rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
# FeaturePlot(object = tmp, features.plot = second.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
# dev.off()

pdf(paste("Marker_expression_Tcells_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("CD8A", "CD4", "IL7R", "CD3G", "CD7", "CD3D", "NKG7", "PTPRC", "NCR1", "NCAM1", "FOXP3")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# marker plot - B cells
pdf(paste("Marker_expression_Bcells_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("CD79A", "CD79B", "BLK", "MS4A1", "CD19", "CXCR4")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# marker plot - Macrophages/DC cells
pdf(paste("Marker_expression_Macrophages_DC_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("CD14", "CD68", "CD163", "CSF1R", "ITGAX", "FCER1A", "CST3")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# marker plot - Endothelial cells
pdf(paste("Marker_expression_Endothelial_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("VWF", "CDH5", "PLVAP")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# Endothelial2
pdf(paste("Marker_expression_Endothelial2_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("VWF", "CDH5", "PLVAP", "SELE", "FLT1", "LYVE1")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# marker plot - CAFs cells
pdf(paste("Marker_expression_CAFs_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("FAP", "THY1", "DCN", "COL1A1", "COL1A2", "COL6A2")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# marker plot - Epithelial cells
pdf(paste("Marker_expression_Epithelial_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("ALCAM", "CD44", "DDR1", "APOC1") # CD166=ALCAM, CD167=DDR1
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

#Epithelial2 cells
pdf(paste("Marker_expression_Epithelial2_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c("UMOD", "SLC12A1", "SPP1", "CA12", "ALDOB", "CALB1", "NAT8", "SLC22A6") # CD166=ALCAM, CD167=DDR1
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# marker plot - Erythrocytes
pdf(paste("Marker_expression_erythrocytes_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c( "HBB", "HBD" ) 
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()


# marker plot - STEM cells
pdf(paste("Marker_expression_stem_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c( "ATXN1", "KIT", "SMO", "ALDH1A1", "UBE2C", "NGFR", "KDM5B" ) # sca-1=ATXN1, c-kit=KIT, cd271=NGFR, JARID1=KDM5B
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

#STEM-2 cells
pdf(paste("Marker_expression_stem2_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c( "NANOG", "POU5F1", "SOX2", "KLF4", "MYC" ) # sca-1=ATXN1, c-kit=KIT, cd271=NGFR, JARID1=KDM5B
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# Other STEM-3 cells
pdf(paste("Marker_expression_stem_other_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c( "ALDH1A1", "POU5F1", "PROM1", "CXCR4" ) # sca-1=ATXN1, c-kit=KIT, cd271=NGFR, JARID1=KDM5B, CD133-PROM1, ALDH1=ALDH1A1, OCT4=POU5F1
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# Driver Genes
pdf(paste("Marker_expression_driverGenes_in_sample_", sample_id, ".pdf", sep=""))
top.markers <- c( "PBRM1", "SETD2", "BAP1", "PTEN", "MET", "KDM5C", "MTOR", "PIK3CA", "TP53" ) 
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()


# new TSNE with cell type
# panc = readRDS( paste("panc_backup_object_in_sample_", sample_id, ".rds", sep="") )

current.cluster.ids <- 0:12
new.cluster.ids <- c("CD4+ T-cells", "NK", "CD8+ T-cells", "Macrophage", "Epithelial(rc, stem)", "DC", "Tumor(Megalin)", "Tumor(AQP1)", "Astrocytes", "Tumor(resistance)", "B-cells", "pro B-cells", "Mesangial")
panc@ident <- plyr::mapvalues(x = panc@ident, from = current.cluster.ids, to = new.cluster.ids)
pdf(paste("TSNE_celltype_in_sample_", sample_id, ".pdf", sep=""))
TSNEPlot(object = panc, do.label = TRUE, pt.size = 1, label.size = 6, no.legend = TRUE, colors.use=c(brewer.pal(9, "Set1"), brewer.pal(4, "Set3")))
dev.off()

# percentage for each cell type

blank_theme <- theme_minimal() +
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

pdf(paste("Cell_type_percentage_in_sample_", sample_id, ".pdf", sep=""))
tmp = as.data.frame(table(panc@ident))
myLegend = paste(tmp$Var1, ": ", tmp$Freq, " (",  round(tmp$Freq/sum(tmp$Freq), digits=3)*100, "%)", sep="")

ggplot(tmp, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c(brewer.pal(9, "Set1"), brewer.pal(4, "Set3")), name=paste("Cell type", sum(tmp$Freq), sep=": "), labels=myLegend) +
  blank_theme +
  theme(legend.position="right")
dev.off()


# save
saveRDS(panc, file = paste("panc_backup_object_in_sample_", sample_id, ".rds", sep=""))



##################################################
##### MONOCLE
##################################################

# look at clusters 4,6,7,9

## create new dataset for monocle
ductal.barcode <- names(panc@ident[panc@ident %in% c(4,6,7,9)])
ductal.raw.data <- panc@raw.data[, ductal.barcode]

fd <- fData(input)
colnames(fd)[2] <- "gene_short_name"
pd <- data.frame(barcode=colnames(ductal.raw.data))
rownames(pd) <- colnames(ductal.raw.data)

HSMM <- newCellDataSet(ductal.raw.data,
                       phenoData = new("AnnotatedDataFrame", data = pd),
                       featureData = new("AnnotatedDataFrame", data = fd),
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

## estimate size factor
HSMM <- estimateSizeFactors(HSMM)

## keep original cluster 
pData(HSMM)$Cluster <- panc@ident[panc@ident %in% c(4,6,7,9)]
pData(HSMM)$Cluster <- droplevels(pData(HSMM)$Cluster)

## consider top DE genes from these clusters of interest
DE_data = panc.markers
DE_data = DE_data[DE_data$cluster==c(4,6,7,9),]
DE_data_subset = DE_data[DE_data$avg_logFC>0.8 & DE_data$p_val_adj<0.00001,]

all.markers.genes.ens = row.names(DE_data_subset)

## genes that define progess
ordering_genes <- all.markers.genes.ens
HSMM <- setOrderingFilter(HSMM, ordering_genes)
sum(fData(HSMM)$use_for_ordering)

## dimension reduction
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')

## ordering the cells in pseudotime
HSMM <- orderCells(HSMM)

## plot cell trajectory
pdf(paste("CellTrajectory_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
plot_cell_trajectory(HSMM, color_by = "Cluster")
dev.off()

pdf(paste("CellTrajectory_in_sample_state_",sample_id, ".pdf", sep=""), width=15, height=9)
plot_cell_trajectory(HSMM, color_by = "State")
dev.off()

###############################
### DE in pseudotime

head(pData(HSMM))

my_pseudotime_de <- differentialGeneTest(HSMM,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 32)

my_pseudotime_de %>% arrange(qval) %>% head()

# save the top 10 genes
my_pseudotime_de %>% arrange(qval) %>% head(n=10L) %>% select(id) -> my_pseudotime_gene
my_pseudotime_gene <- my_pseudotime_gene$id

# PLOT top genes change in pseudotime
pdf( paste( 'top_differential_genes_in_pseudotime_top10_by_state_', sample_id, '.pdf', sep="")  )
plot_genes_in_pseudotime(HSMM[my_pseudotime_gene,])
dev.off()

pdf( paste( 'top_differential_genes_in_pseudotime_top10_by_cluster_', sample_id, '.pdf', sep="")  )
plot_genes_in_pseudotime(HSMM[my_pseudotime_gene,], color_by = "Cluster")
dev.off()

###########################################
######## DEEP ANALYSIS
###########################################

panc = readRDS( paste("panc_backup_object_in_sample_", sample_id, ".rds", sep="") )

# list of cell barcodes in the "Tumor(AQP1)" cluster
cells_to_take = names(panc@ident)[panc@ident=="Tumor(AQP1)"]

# make a new object using the list of cells above
# new_panc = FilterCells(object = panc, cells.use = cells_to_take )
# new_panc = FetchData( panc, cells.use=cells_to_take, use.scaled = T )

new_panc = SubsetData(object = panc, ident.use = "Tumor(AQP1)", subset.raw = T, do.clean = T)

new_panc = run_scale_pca_cluster( new_panc, "Tumor_AQP1" )

epithelial_markers = c("UMOD", "SLC12A1", "SPP1", "CA12", "ALDOB", "CALB1", "PDZK1IP1", "NAT8", "SLC22A6", "SLC22A8", 
                       "SLC34A1", "SLC12A3", "SCNN1B", "CLCN5", "CLDN16", "GPX3", "DEFB1", "KCNJ1", "KNG1", "SCNN1A", 
                       "SLC22A2", "PAX8", "SLC23A3", "KCNJ15", "MT1G", "SLC12A6", "BHMT", "ALDOB", "SLC13A1")

endothelial_markers = c("SELE", "FLT1", "LYVE1", "VWF", "MCAM", "KDR", "CDH5", "ARHGDIB", "A2M", "PTPRB") # FLK1=KDR
tumor_markers = c("NPHS2", "LRP2", "AQP1", "AQP2", "UMOD")
new_pc_genes = c(epithelial_markers, endothelial_markers, tumor_markers)

new_pc_genes_ens <- fData(input)[match(new_pc_genes, fData(input)$symbol), 1]

new_panc <- RunPCA(new_panc, pc.genes = new_pc_genes_ens, do.print = F)

new_panc <- FindClusters(new_panc, reduction.type = "pca", dims.use = 1:10, resolution = 0.8, print.output = 0, save.SNN = TRUE, force.recalc=T)

new_panc <- RunTSNE(new_panc, dims.use = 1:10, do.fast = TRUE)

name = "Tumor_AQPT2"

# save
saveRDS(new_panc, file = paste("new_panc_backup_object_in_sample_", name, "_", sample_id, ".rds", sep=""))

pdf( paste( 'TSNE_', name, "_", sample_id, '.pdf', sep="")  )
TSNEPlot(new_panc, pt.size=2, label.size = 6)
dev.off()

# Endothelial2
pdf(paste("Marker_expression_Endothelial2_", name, "_", sample_id, ".pdf", sep=""))
top.markers <- c("VWF", "CDH5", "PLVAP", "SELE", "FLT1", "LYVE1")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- new_panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

#Epithelial2 cells
pdf(paste("Marker_expression_Epithelial2_", name, "_", sample_id, ".pdf", sep=""))
top.markers <- c("UMOD", "SLC12A1", "SPP1", "CA12", "ALDOB", "CALB1", "NAT8", "SLC22A6") # CD166=ALCAM, CD167=DDR1
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- new_panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# STEM
pdf(paste("Marker_expression_stem2_in_sample_", name, "_", sample_id, ".pdf", sep=""))
top.markers <- c( "NANOG", "POU5F1", "SOX2", "KLF4", "MYC" ) # sca-1=ATXN1, c-kit=KIT, cd271=NGFR, JARID1=KDM5B
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- new_panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

pdf(paste("Marker_expression_stem3_other_in_sample_", name, "_", sample_id, ".pdf", sep=""))
top.markers <- c( "ALDH1A1", "POU5F1", "PROM1", "CXCR4" ) # sca-1=ATXN1, c-kit=KIT, cd271=NGFR, JARID1=KDM5B, CD133-PROM1, ALDH1=ALDH1A1, OCT4=POU5F1
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- new_panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()

# tumor marker plot
pdf(paste("Tumor_marker_expression_", name, "_", sample_id, ".pdf", sep=""))
top.markers <- c("NPHS2", "LRP2", "AQP1", "AQP2", "UMOD")
markers_ens <- fData(input)[match(top.markers, fData(input)$symbol), 1]
tmp <- new_panc
tmp@data <- tmp@data[rownames(tmp@data) %in% markers_ens, ]
rownames(tmp@data) <- fData(input)[match(rownames(tmp@data), fData(input)$id), 2]
FeaturePlot(object = tmp, features.plot = top.markers, cols.use = c("grey", "red"), reduction.use = "tsne")
dev.off()




