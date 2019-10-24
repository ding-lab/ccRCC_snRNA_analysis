# Yige Wu @WashU Sep 2019
## for integrating two snRNA datasets for sample CPT0086820004 and CPT0075130004 (from cellranger output with premrna reference)

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# set parameters ----------------------------------------------------------
version_tmp <- 3
sample_id <- "CPT0086820004_CPT007513000"

###########################################
######## Dataset preprocessing
###########################################

# Input seurat objects -----------------------------------------------------
## v1, remove the log normalizatin step and find variable feature step, still the same
## v2 changes the orders of UMAP and find neighbors, changes the dimension number, parameters starting from ScaleData
## v3, add log normalization step back to v2
CPT0086820004_obj <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/CPT0086820004/CPT0086820004/Seurat_outs/pf5000_low.nFeature200_high.nFeature10000_low.nCount5000_high.nCount10000_high.percent.mito0.1/single_cell_study_processed.rds")
CPT0086820004_obj$orig.ident  <- "CPT0086820004"
DefaultAssay(CPT0086820004_obj) <- "RNA" #only have 3000 features
CPT0086820004_obj@assays$SCT <- NULL
CPT0086820004_obj@graphs <- list()
CPT0086820004_obj@neighbors <- list()
CPT0086820004_obj@reductions <- list()
for (col_name_tmp in c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.5", "seurat_clusters")) {
  if (col_name_tmp %in% names(CPT0086820004_obj@meta.data)) {
    CPT0086820004_obj@meta.data[[col_name_tmp]] <- NULL
  }
}

CPT0075130004_obj <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/CPT0075130004/CPT0075130004/Seurat_outs/pf5000_low.nFeature200_high.nFeature10000_low.nCount5000_high.nCount10000_high.percent.mito0.1/single_cell_study_processed.rds")
CPT0075130004_obj$orig.ident  <- "CPT0075130004"
DefaultAssay(CPT0075130004_obj) <- "RNA" #only have 3000 features
CPT0075130004_obj@assays$SCT <- NULL
CPT0075130004_obj@graphs <- list()
CPT0075130004_obj@neighbors <- list()
CPT0075130004_obj@reductions <- list()
for (col_name_tmp in c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.5", "seurat_clusters")) {
  if (col_name_tmp %in% names(CPT0075130004_obj@meta.data)) {
    CPT0075130004_obj@meta.data[[col_name_tmp]] <- NULL
  }
}

# Merge seurat objects -----------------------------------------------------
renal_merged <- merge(x = CPT0086820004_obj, y = CPT0075130004_obj, 
                      add.cell.ids = c('CPT0086820004','CPT0075130004'),
                      project = 'Kidney')

#  split the combined object into a list, with each dataset as an element ----------------------------------------
renal.list <- SplitObject(renal_merged, split.by = "orig.ident")
renal.list <- lapply(X = renal.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


###########################################
######## Integration of datasets
###########################################

# create reference anchors  ----------------------------------------------------
## Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input. 
## Here, we integrate three of the objects into a reference (we will use the fourth later in this vignette)
reference.list <- renal.list[c("CPT0086820004", "CPT0075130004")]
renal.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)
## We use all default parameters here for identifying anchors, including the ‘dimensionality’ of the dataset (30; feel free to try varying this parameter over a broad range, for example between 10 and 50).


# Integrate Data -----------------------------------------------------------
## We then pass these anchors to the IntegrateData function, which returns a Seurat object.
renal.integrated <- IntegrateData(anchorset = renal.anchors, dims = 1:20)
## The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
## After running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix. Note that the original (uncorrected values) are still stored in the object in the “RNA” assay, so you can switch back and forth.

## We can then use this new integrated matrix for downstream analysis and visualization. Here we scale the integrated data, run PCA, and visualize the results with UMAP. 


# switch to integrated assay. ---------------------------------------------
DefaultAssay(renal.integrated) <- "integrated" #only have 3000 features

# Run the standard workflow for visualization and clustering ------------
# The variable features of this assay are automatically set during IntegrateData
renal.integrated <- ScaleData(renal.integrated, verbose = F) 
renal.integrated <- RunPCA(renal.integrated, npcs = 30, verbose = FALSE)
renal.integrated <- RunUMAP(renal.integrated, reduction = "pca", dims = 1:20)
renal.integrated <- FindNeighbors(renal.integrated, reduction = "pca", dims = 1:20)
renal.integrated <- FindClusters(renal.integrated, resolution = 0.5)
saveRDS(object = renal.integrated, file = paste0(makeOutDir(), sample_id, "_renal_integrated.20190909.v1.RDS"))

# Plot the dimension reduction grouped by sample and grouped by cluster -------------
library(cowplot)
## make sure the grouping variable is in the meta data
renal.integrated@meta.data %>%
  head()

file2write <- paste(makeOutDir(), "DimPlot_raw_cell_clusters_in_sample_", sample_id, ".20190904.", "v", version_tmp, ".pdf", sep="")
pdf(file = file2write, width = 15, height = 6)
p1 <- DimPlot(renal.integrated, reduction = "umap", group.by = "orig.ident", label = T)
p2 <- DimPlot(renal.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)
plot_grid(p1, p2)
dev.off()

file2write <- paste(makeOutDir(), "DimPlot_raw_cell_clusters_in_sample_", sample_id, ".splitted.20190904.", "v", version_tmp, ".pdf", sep="")
pdf(file = file2write, width = 12, height = 6)
DimPlot(renal.integrated, reduction = "umap", split.by = "orig.ident", label = T)
dev.off()

###########################################
######## Differential expression
###########################################

# find DEG ----------------------------------------------------------------
renal.integrated <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/scRNA/integrate_CPT0086820004_CPT0075130004_20190904_v2/CPT0086820004_CPT007513000_renal_integrated.20190909.v1.RDS")
DefaultAssay(renal.integrated) <- "RNA"

renal.markers <- FindAllMarkers(object = renal.integrated, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
renal.markers %>%
  colnames()
renal.markers <- renal.markers[, c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")]
write.table(renal.markers, file = paste0(makeOutDir(), "Renal.Markers.Pos.txt"), quote = F, sep = "\t", row.names = F)


# filter DEG by manual curated markers ------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Kidney_Markers/RCC_literature_review_and_marker_gene_table - Gene_to_Cell_Type_Table.20190911.v2.tsv")
for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr)) {
  marker_genes_tmp <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  renal.markers[, cell_type_tmp] <- as.numeric(renal.markers$gene %in% marker_genes_tmp)
}
write.table(renal.markers, file = paste0(makeOutDir(), "Renal.Markers.Pos.CellTypeMarkerAnnotated.txt"), quote = F, sep = "\t", row.names = F)

renal.markers.filtered <- renal.markers[rowSums(renal.markers[,unique(gene2cellType_tab$Cell_Type_Abbr)]) > 0,]

cell_types2print <- unique(gene2cellType_tab$Cell_Type_Abbr)[colSums(renal.markers[,unique(gene2cellType_tab$Cell_Type_Abbr)]) > 0]
renal.markers.filtered <- renal.markers.filtered[, c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene", cell_types2print)]

write.table(renal.markers.filtered, file = paste0(makeOutDir(), "renal.markers.pos.filtered.txt"), quote = F, sep = "\t", row.names = F)

renal.markers %>% head()
DEdata = renal.markers %>% group_by(cluster) %>% top_n(100, avg_logFC)
write.table(DEdata, file = paste0(makeOutDir(), "Top100_DE_only.pos.txt"), quote = F, sep = "\t", row.names = F)

renal.markers.bycluster <- renal.markers %>%
  filter(cluster == 0)

## Confident cluster 0: VCAM1+ ccRCC #1
### Defining markers:
#### ccRCC: CA9, LRP2
#### Proximal convoluted: SLC17A3
#### Proximal straight: SLC16A9
### Confounding markers: 
#### Vascular: VCAM1
#### CapMesenchyme: PAX2
#### UretericBud: HNF1B

## Cluster 1: NK
### Defining markers:
#### NK: KLRD1
### Confounding markers: 
#### CapMesenchyme: PAX2
#### UretericBud: TFCP2L1

## Cluster 2: Tumor associated Macrophage
### Defining Postive markers:
#### TAM: FCGR3A; CD163; CSF1R; CSF3R; CD86; MRC1;
### Defining Negative markers: CD3-
### Confounding markers: 
#### Vascular: PECAM1
#### NK: TYROBP (not so specific),
#### UretericBud: HNF1B

## Cluster 3: VCAM+ Mast Cell
### Defining markers:
#### Mast Cell: ENPP3
### Confounding markers: 
#### Vascular: VCAM1

## Uncertain Cluster 4: 
### Defining markers:
#### Convoluted Proximal Tubule: SLC17A3
### Confounding markers: 
#### NK: KLRD1 (check literature)

## Uncertain Cluster 5: 
### Defining markers:
### Confounding markers: 
#### Vascular: PECAM1
#### CD4Tcell: CD4
#### Bcell: HLA-DRA

## Cluster 6: MET+ Proximal Tubule
### Defining markers:
#### Proximal Tubule: SLC13A3
#### pRCC: MET
### Confounding markers: 
#### Mast Cell: ENPP3
#### NK: NCAM1

## cluster 7: Proliferating CD8+ cell
### Defining markers:
#### Proximal Tubule: SLC13A3
#### pRCC: MET
### Confounding markers: 
#### Mast Cell: ENPP3
#### NK: NCAM1

## Cluster 8: MesangialCell/Myofibroblast
### Defining markers:
#### MesangialCell: PDGFRB, ACTA2, NDUFA4L2
### Confounding Markers
#### Endothelial: VWF
#### AscendingVasaRecta: PLVAP, PTPRB

## Cluster 9: Endothelial & Podocyte
### Defining markers:
#### AscendingVasaRecta: PLVAP, KDR, PTPRB, PECAM1
#### Podocyte: PODXL
### Confounding Markers
#### DC: NRP1

## Uncertain cluster 10:
### Confounding Markers
#### Podocyte: PTPRO
#### CD4Tcell: CD4
#### Neutrophil: MNDA
#### Bcell: HLA-DRA
#### DC: NRP1
#### Vascular: PECAM1
#### Proliferating: MKI67, TOP2A

## Confident Cluster 11: proliferating ccRCC #2
### Defining markers:
#### ccRCC: CA9, LRP2
#### Proximal convoluted: SLC17A3
#### Proliferating: MKI67, TOP2A

## Uncertain Cluster 12:
#### Podocyte: PTPRO
#### CD4Tcell: CD4
#### Neutrophil: MNDA
#### Bcell: HLA-DRA, CD40
#### DC: IRF8
#### NK: NCAM1

## Cluster 13: NK
### Defining markers:
#### NK: FCGR3A, TYROBP
### Confounding markers: 
#### Vascular: PECAM1
#### Tcell: IL32
#### Neutrophil: MNDA

# DEG between tumor clusters ----------------------------------------------


###########################################
######## Plotting marker expression
###########################################
object2plot <- renal.integrated
DefaultAssay(object2plot) <- "RNA"

outpath <- makeOutDir()
marker_exp_out_path <- paste0(outpath, "Marker_Expression/")
dir.create(marker_exp_out_path)

plot_parameters <- list()

plot_parameters[["Podocyte"]] <- list()
plot_parameters[["Podocyte"]][["top.markers"]] <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == "Podocyte"]
plot_parameters[["Podocyte"]][["FeaturePlot"]] <- list() 
plot_parameters[["Podocyte"]][["FeaturePlot"]][["width"]] <- 12
plot_parameters[["Podocyte"]][["FeaturePlot"]][["height"]] <- 16
plot_parameters[["Podocyte"]][["VlnPlot"]] <- list() 
plot_parameters[["Podocyte"]][["VlnPlot"]][["width"]] <- 6
plot_parameters[["Podocyte"]][["VlnPlot"]][["height"]] <- 16


plot_parameters[["TAM"]] <- list()
plot_parameters[["TAM"]][["top.markers"]] <- c("CD163", "CSF1R", "CSF3R", "CD68", "CD86", "MRC1")
plot_parameters[["TAM"]][["FeaturePlot"]] <- list() 
plot_parameters[["TAM"]][["FeaturePlot"]][["width"]] <- 14
plot_parameters[["TAM"]][["FeaturePlot"]][["height"]] <- 16
plot_parameters[["TAM"]][["VlnPlot"]] <- list() 
plot_parameters[["TAM"]][["VlnPlot"]][["width"]] <- 6
plot_parameters[["TAM"]][["VlnPlot"]][["height"]] <- 16

plot_parameters[["ccRCC"]] <- list()
plot_parameters[["ccRCC"]][["top.markers"]] <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == "ccRCC"]
plot_parameters[["ccRCC"]][["FeaturePlot"]] <- list() 
plot_parameters[["ccRCC"]][["FeaturePlot"]][["width"]] <- 12
plot_parameters[["ccRCC"]][["FeaturePlot"]][["height"]] <- 16
plot_parameters[["ccRCC"]][["VlnPlot"]] <- list() 
plot_parameters[["ccRCC"]][["VlnPlot"]][["width"]] <- 6
plot_parameters[["ccRCC"]][["VlnPlot"]][["height"]] <- 16

plot_parameters[["Proximal_Tubule"]] <- list()
plot_parameters[["Proximal_Tubule"]][["top.markers"]] <- c("SLC22A8", "SLC17A3" , "SLC22A7", "SLC16A9", "SLC7A13", "SLC34A1", "SLC13A3")
plot_parameters[["Proximal_Tubule"]][["FeaturePlot"]] <- list() 
plot_parameters[["Proximal_Tubule"]][["FeaturePlot"]][["width"]] <- 14
plot_parameters[["Proximal_Tubule"]][["FeaturePlot"]][["height"]] <- 16
plot_parameters[["Proximal_Tubule"]][["VlnPlot"]] <- list() 
plot_parameters[["Proximal_Tubule"]][["VlnPlot"]][["width"]] <- 6
plot_parameters[["Proximal_Tubule"]][["VlnPlot"]][["height"]] <- 16

plot_parameters[["RCC_Oncogene"]] <- list()
plot_parameters[["RCC_Oncogene"]][["top.markers"]] <- c("HIF1A", "EPAS1" , "MTOR", "MET", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3")
plot_parameters[["RCC_Oncogene"]][["FeaturePlot"]] <- list() 
plot_parameters[["RCC_Oncogene"]][["FeaturePlot"]][["width"]] <- 16
plot_parameters[["RCC_Oncogene"]][["FeaturePlot"]][["height"]] <- 16
plot_parameters[["RCC_Oncogene"]][["VlnPlot"]] <- list() 
plot_parameters[["RCC_Oncogene"]][["VlnPlot"]][["width"]] <- 6
plot_parameters[["RCC_Oncogene"]][["VlnPlot"]][["height"]] <- 20

plot_parameters[["RCC_TSG"]] <- list()
plot_parameters[["RCC_TSG"]][["top.markers"]] <- c("VHL", "PBRM1" , "BAP1", "SETD2", "PTEN", "TP53")
plot_parameters[["RCC_TSG"]][["FeaturePlot"]] <- list() 
plot_parameters[["RCC_TSG"]][["FeaturePlot"]][["width"]] <- 14
plot_parameters[["RCC_TSG"]][["FeaturePlot"]][["height"]] <- 16
plot_parameters[["RCC_TSG"]][["VlnPlot"]] <- list() 
plot_parameters[["RCC_TSG"]][["VlnPlot"]][["width"]] <- 6
plot_parameters[["RCC_TSG"]][["VlnPlot"]][["height"]] <- 16


# plotting ----------------------------------------------------------------
for (cell_type_tmp in names(plot_parameters)) {
  file2write <- paste0(marker_exp_out_path, cell_type_tmp, "_marker_expression_in_sample_", sample_id, "_FeaturePlot.pdf")
  pdf(file2write, 
      width = plot_parameters[[cell_type_tmp]][["FeaturePlot"]][["width"]],
      height = plot_parameters[[cell_type_tmp]][["FeaturePlot"]][["height"]])
  p <- FeaturePlot(object = object2plot, features = plot_parameters[[cell_type_tmp]][["top.markers"]], 
                   cols = c("grey", "red"), reduction = "umap", label = T)
  print(p)
  dev.off()
  
  pdf(paste(marker_exp_out_path, cell_type_tmp, "_marker_expression_in_sample_", sample_id, "_VlnPlot.pdf", sep=""), 
      width = plot_parameters[[cell_type_tmp]][["VlnPlot"]][["width"]],
      height = plot_parameters[[cell_type_tmp]][["VlnPlot"]][["height"]])
  p <- VlnPlot(object2plot, features = plot_parameters[[cell_type_tmp]][["top.markers"]],
               ncol = 1, pt.size = 0.5)
  print(p)
  dev.off()
}


# plot immnue cell markers ----------------------------------------------------------
outpath <- makeOutDir()
object2plot <- renal.integrated
marker_exp_out_path <- paste0(outpath, "Immune_Marker_Expression/")
dir.create(marker_exp_out_path)

Plasma_markers <- data.frame(gene_symbol = c("SDC1", "IGHG1", "IGHG3", "IGHG4"), cell_type = "Plasma")
B_markers <- data.frame(gene_symbol = c("CD19", "MS4A1", "CD79A", "CD79B"), cell_type = "B")
Monocyte_markers <- data.frame(gene_symbol = c("LYZ", "CD14", "S100A8"), cell_type = "Monocyte")
Macrophage_markers <- data.frame(gene_symbol = c("FCGR3A", "MS4A7", "IFITM3"), cell_type = "Macrophage")
DC_markers <- data.frame(gene_symbol = c("FCER1A"), cell_type = "DC")
CD8T_markers <- data.frame(gene_symbol = c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B"), cell_type = "CD8+T")
CD4T_markers <- data.frame(gene_symbol = c("CD3D", "CD3E", "CD3G", "IL7R", "LDHB", "NOSIP", "CD4"), cell_type = "CD4+T")
NK_markers <- data.frame(gene_symbol = c("GNLY", "TYROBP", "HOPX", "FCGR3A"), cell_type = "NK")
Erythrocytes_markers <- data.frame(gene_symbol = c("HBD", "GYPA", "HBA1", "HBA2", "CA1"), cell_type = "Erythrocyte")

## plot CD4Tcell Markers
pdf(paste(marker_exp_out_path,"CD4Tcell_Marker_expression_in_sample_", sample_id, ".pdf", sep=""),height=12,width=15)
top.markers <- as.vector(CD4T_markers$gene_symbol)
FeaturePlot(object = object2plot, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

## plot CD8Tcell Markers
pdf(paste(marker_exp_out_path,"CD8Tcell_Marker_expression_in_sample_", sample_id, ".pdf", sep=""),height=15,width=12)
top.markers <- as.vector(CD8T_markers$gene_symbol)
FeaturePlot(object = object2plot, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

## plot Macrophage Markers
pdf(paste(marker_exp_out_path,"Macrophage_Marker_expression_in_sample_", sample_id, ".pdf", sep=""),height=12,width=12)
top.markers <- as.vector(Monocyte_markers$gene_symbol)
FeaturePlot(object = object2plot, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

## plot Monocyte Markers
pdf(paste(marker_exp_out_path,"Monocyte_Marker_expression_in_sample_", sample_id, ".pdf", sep=""),height=6,width=8)
top.markers <- as.vector(Macrophage_markers$gene_symbol)
FeaturePlot(object = object2plot, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

## plot NK Markers
pdf(paste(marker_exp_out_path,"NK_Marker_expression_in_sample_", sample_id, ".pdf", sep=""),height=6,width=8)
top.markers <- as.vector(NK_markers$gene_symbol)
FeaturePlot(object = object2plot, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

## plot Plasma Markers
pdf(paste(marker_exp_out_path,"Plasma_Marker_expression_in_sample_", sample_id, ".pdf", sep=""),height=12,width=12)
top.markers <- as.vector(Plasma_markers$gene_symbol)
FeaturePlot(object = object2plot, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

## plot Erythrocytes Markers
pdf(paste(marker_exp_out_path,"Erythrocytes_Marker_expression_in_sample_", sample_id, ".pdf", sep=""),height=12,width=12)
top.markers <- as.vector(Erythrocytes_markers$gene_symbol)
FeaturePlot(object = object2plot, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

## plot BCell Markers
pdf(paste(marker_exp_out_path,"BCells_Marker_expression_in_sample_", sample_id, ".pdf", sep=""),height=12,width=12)
top.markers <- as.vector(B_markers$gene_symbol)
FeaturePlot(object = object2plot, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

## plot DC Markers
pdf(paste(marker_exp_out_path,"DC_Marker_expression_in_sample_", sample_id, ".pdf", sep=""),height=12,width=12)
top.markers <- as.vector(DC_markers$gene_symbol)
FeaturePlot(object = object2plot, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

# plot Endothelial cell markers ----------------------------------------------------------
marker_exp_out_path <- paste0(outpath, "Endothelial_Marker_Expression/")
dir.create(marker_exp_out_path)

## plot Endomthelial markers
pdf(paste(marker_exp_out_path,"Endothelial_marker_expression_in_sample_", sample_id, ".pdf", sep=""),height=14,width=12)
top.markers <- c("VWF", "PECAM1", "FLT4", "FLT1", "FLT3", "KDR")
FeaturePlot(object = object2plot, features = top.markers, cols = c("grey", "red"), reduction = "umap")
dev.off()

###########################################
######## Label Cells
###########################################

# Label Cells -  --------
n_clusters <- CPT0086820004_obj@meta.data$seurat_clusters %>% nlevels
current.cluster.ids <- 0:(n_clusters-1)
new.cluster.ids<- c('CD8Tcell','CD4Tcell','Macrophage','NK(1)',
                    'Epithelial(rc, stem)','NK(2)','DC','NK(3)','Tumor(Megalin)',
                    'Treg','NK(4)','Endothelial','Astrocyte','Tumor(Resistance)',
                    'Bcell','Pro-Bcell','Unknown')
renal_a1@meta.data$cell_type <- plyr::mapvalues(as.character(Idents(renal_a1)), from= current.cluster.ids,to = new.cluster.ids)
Idents(renal_a1) <- renal_a1@meta.data$cell_type
table(renal_a1@meta.data$cell_type)


# trasnfer the labels -----------------------------------------------------
###### label transfer: A1 to integrated
renal.anchors <- FindTransferAnchors(reference = CPT0086820004_obj, query = renal.integrated, 
                                     dims = 1:30)
predictions <- TransferData(anchorset = renal.anchors, refdata = CPT0086820004_obj$cell_type, 
                            dims = 1:30)
renal.integrated@meta.data$cell_type <- predictions$predicted.id

Idents(renal.integrated) <- renal.integrated@meta.data$cell_type

table(renal.integrated@meta.data$cell_type)
DimPlot(renal.integrated, label = TRUE) + NoLegend()
