library(cellrangerRkit)
library(Seurat)
library(dplyr)

# input
sample_id <- "TWAN-RCC-B2"
cell_path <- "/diskmnt/Datasets/Kidney/scRNA"
# cell_path = "~/DingLab/SingleCell/Data/Kidney_Data/"
input <- get_matrix_from_h5(paste(cell_path, "/", sample_id, "/outs/raw_gene_bc_matrices_h5.h5", sep=""))

# pre-filter
bc_sums <- colSums(input)
bg_id <- names(bc_sums[bc_sums < 300])

filter_func <- function(gbm, id_to_remove)
{
  new_matrix <- exprs(gbm)[,!colnames(exprs(gbm)) %in% id_to_remove]
  pd <- data.frame(barcode = colnames(new_matrix))
  rownames(pd) <- colnames(new_matrix)
  return(newGeneBCMatrix(new_matrix, fd = fData(gbm), pd = pd, template = gbm))
}

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
panc <- FilterCells(object = panc, subset.names = c("nGene", "nUMI", "percent.mito"), low.thresholds = c(200, 1000, -Inf), high.thresholds = c(7500, 50000, 0.10))
nrow(FetchData(panc))
# TWAN-RCC-B2: 396
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
DEdata = panc.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
write.table(DEdata, file = "./DEdata.txt", quote = F, sep = "\t")
# cluster0 marker gene:IL7R, CD8A, CCD7(T cells- consistts of cd4+, cd8+ and NK)
# cluster1 marker gene:ITGAX (DC, has some MPP signature TIMP1)
# cluster2 marker gene:APOC1 (Macrophages)
# cluster3 marker gene:CRYAB(Endothelial + bit tumor from LRP2)


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

# new TSNE with cell type

current.cluster.ids <- 0:3
new.cluster.ids <- c("T-cells", "DC", "Macrophage", "Epithelial(rc, stem)")
panc@ident <- plyr::mapvalues(x = panc@ident, from = current.cluster.ids, to = new.cluster.ids)
pdf(paste("TSNE_celltype_in_sample_", sample_id, ".pdf", sep=""))
TSNEPlot(object = panc, do.label = TRUE, pt.size = 1, label.size = 6, no.legend = TRUE, colors.use=brewer.pal(4, "Set1"))
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
  scale_fill_manual(values=brewer.pal(4, "Set1"), name=paste("Cell type", sum(tmp$Freq), sep=": "), labels=myLegend) +
  blank_theme +
  theme(legend.position="right")
dev.off()


# save
saveRDS(panc, file = paste("panc_backup_object_in_sample_", sample_id, ".rds", sep=""))
