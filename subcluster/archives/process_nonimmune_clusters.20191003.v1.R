# Yige Wu @WashU Sep 2019
## for isolating the non-immune cell clusters and re-do clustering

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
sample_ids <- c("CPT0075140002", "CPT0001260013", "CPT0086350004")

# Input the seurat object -------------------------------
renal.int.object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/integration/integrate_seurat_objects/20191002.v1/renal_integrated20191002.v1.RDS")

# input cell type assignment table ----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/Cell_Type_Assignment/RCC_snRNA_Downstream_Processing - Cluster2Cell_Type_Assignment.20191002.v1.tsv", data.table = F)
nonimmune_clusters <- cluster2celltype_tab %>%
  filter(Is_Immune == "No") %>%
  select(Cluster)

nonimmune_clusters <- nonimmune_clusters$Cluster

# input the gene2cell type annotation file --------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191003.v1.tsv")

###########################################
######## Dataset preprocessing
###########################################

# subset object by non-immune clusters ------------------------------------
## subset usage: https://satijalab.org/seurat/v3.0/interaction_vignette.html
renal.int.object <- subset(renal.int.object, idents = nonimmune_clusters)
renal.int.object@meta.data %>% head()


# create list of seurat objects by sample ---------------------------------
renal.list <- list()
# Input seurat objects -----------------------------------------------------
for (sample_id_tmp in sample_ids) {
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

###########################################
######## Integration of datasets
###########################################

# create reference anchors  ----------------------------------------------------
## Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input. 
## Here, we integrate three of the objects into a reference (we will use the fourth later in this vignette)
renal.anchors <- FindIntegrationAnchors(object.list = renal.list, dims = 1:20)
## We use all default parameters here for identifying anchors, including the ‘dimensionality’ of the dataset (30; feel free to try varying this parameter over a broad range, for example between 10 and 50).


# Integrate Data -----------------------------------------------------------
## We then pass these anchors to the IntegrateData function, which returns a Seurat object.
renal.int.nonimmune.obj <- IntegrateData(anchorset = renal.anchors, dims = 1:20)
## The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
## After running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix. Note that the original (uncorrected values) are still stored in the object in the “RNA” assay, so you can switch back and forth.

## We can then use this new integrated matrix for downstream analysis and visualization. Here we scale the integrated data, run PCA, and visualize the results with UMAP. 


# switch to integrated assay. ---------------------------------------------
DefaultAssay(renal.int.nonimmune.obj) <- "integrated" #only have 3000 features

# Run the standard workflow for visualization and clustering ------------
# The variable features of this assay are automatically set during IntegrateData
renal.int.nonimmune.obj <- ScaleData(renal.int.nonimmune.obj, verbose = F) 
renal.int.nonimmune.obj <- RunPCA(renal.int.nonimmune.obj, npcs = 30, verbose = FALSE)
renal.int.nonimmune.obj <- RunUMAP(renal.int.nonimmune.obj, reduction = "pca", dims = 1:20)
renal.int.nonimmune.obj <- FindNeighbors(renal.int.nonimmune.obj, reduction = "pca", dims = 1:20)
renal.int.nonimmune.obj <- FindClusters(renal.int.nonimmune.obj, resolution = 0.5)
saveRDS(object = renal.int.nonimmune.obj, file = paste0(dir_out, "NonImmune_Integrated.", run_id, ".RDS"))

###########################################
######## Plotting dimension
###########################################

# plot dimention reduction ---------------------------------------------
file2write <- paste(dir_out, "DimPlot_NonImmune_Clusters.", run_id, ".pdf", sep="")
pdf(file = file2write, width = 8, height = 6)
p2 <- DimPlot(renal.int.nonimmune.obj, reduction = "umap", group.by = "seurat_clusters", label = T)
print(p2)
dev.off()

file2write <- paste(dir_out, "DimPlot_NonImmune_Clusters_SampleLabeled.", run_id, ".pdf", sep="")
pdf(file = file2write, width = 15, height = 6)
p1 <- DimPlot(renal.int.nonimmune.obj, reduction = "umap", group.by = "orig.ident", label = F)
p2 <- DimPlot(renal.int.nonimmune.obj, reduction = "umap", group.by = "seurat_clusters", label = T)
plot_grid(p1, p2)
dev.off()

file2write <- paste(dir_out, "DimPlot_NonImmune_Clusters_SampleSplitted.", run_id, ".pdf", sep="")
pdf(file = file2write, width = 16, height = 6)
DimPlot(renal.int.nonimmune.obj, reduction = "umap", split.by = "orig.ident", label = T)
dev.off()

###########################################
######## Plotting marker expression
###########################################
object2plot <- renal.int.nonimmune.obj
DefaultAssay(object2plot) <- "RNA"

dir_marker_exp <- paste0(dir_out, "Cell_Type_Marker_Expression/")
dir.create(dir_marker_exp)

for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr[gene2cellType_tab$Cell_Type_Abbr != "Other"])) {
  markers2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  file2write <- paste0(dir_marker_exp, cell_type_tmp, "_Marker_Expression.", run_id, ".pdf")
  pdf(file2write, 
      width = 12,
      height = 16)
  p <- FeaturePlot(object = object2plot, features = markers2plot, 
                   cols = c("grey", "red"), reduction = "umap", label = T)
  print(p)
  dev.off()
}

dir_marker_exp <- paste0(dir_out, "Other_Marker_Expression/")
dir.create(dir_marker_exp)

for (biomarker_type_tmp in unique(gene2cellType_tab$Biomarker_Type[gene2cellType_tab$Cell_Type_Abbr == "Other"])) {
  markers2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Biomarker_Type == biomarker_type_tmp]
  file2write <- paste0(dir_marker_exp, biomarker_type_tmp, "_Marker_Expression.", run_id, ".pdf")
  pdf(file2write, 
      width = 12,
      height = 16)
  p <- FeaturePlot(object = object2plot, features = markers2plot, 
                   cols = c("grey", "red"), reduction = "umap", label = T)
  print(p)
  dev.off()
}

dir_markers_exp <- paste0(dir_out, "Marker_Expression_by_gene/")
dir.create(dir_markers_exp)

for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr[gene2cellType_tab$Cell_Type_Abbr != "Other"])) {
  markers2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  markers2plot <- markers2plot[markers2plot %in% rownames(object2plot@assays$RNA@counts)]
  
  if (length(markers2plot) > 0) {
    dir_cell_type_markers_exp <- paste0(dir_markers_exp, cell_type_tmp, "/")
    dir.create(dir_cell_type_markers_exp)
    
    for (marker_tmp in markers2plot) {
      file2write <- paste0(dir_cell_type_markers_exp, marker_tmp, "_Marker_Expression.", run_id, ".pdf")
      
      pdf(file2write, 
          width = 7,
          height = 6)
      p <- FeaturePlot(object = object2plot, features = marker_tmp, 
                       cols = c("grey", "red"), reduction = "umap", label = T)
      print(p)
      dev.off()
    }
  }
}

for (biomarker_type_tmp in unique(gene2cellType_tab$Biomarker_Type[gene2cellType_tab$Cell_Type_Abbr == "Other"])) {
  markers2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Biomarker_Type == biomarker_type_tmp]
  markers2plot <- markers2plot[markers2plot %in% rownames(object2plot@assays$RNA@counts)]
  
  if (length(markers2plot) > 0) {
    dir_cell_type_markers_exp <- paste0(dir_markers_exp, "Other_", biomarker_type_tmp, "/")
    dir.create(dir_cell_type_markers_exp)
    
    for (marker_tmp in markers2plot) {
      file2write <- paste0(dir_cell_type_markers_exp, marker_tmp, "_Marker_Expression.", run_id, ".pdf")
      
      pdf(file2write, 
          width = 7,
          height = 6)
      p <- FeaturePlot(object = object2plot, features = marker_tmp, 
                       cols = c("grey", "red"), reduction = "umap", label = T)
      print(p)
      dev.off()
    }
  }
}
rm(object2plot)

###########################################
######## Differential expression
###########################################

# find DEG ----------------------------------------------------------------
DefaultAssay(renal.int.nonimmune.obj) <- "RNA"

renal.markers <- FindAllMarkers(object = renal.int.nonimmune.obj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
renal.markers %>%
  colnames()
renal.markers <- renal.markers[, c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
write.table(renal.markers, file = paste0(dir_out, "NonImmune.DEGs.Pos.txt"), quote = F, sep = "\t", row.names = F)

# filter DEG by manual curated markers ------------------------------------

for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr[gene2cellType_tab$Cell_Type_Abbr != "Other"])) {
  marker_genes_tmp <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  renal.markers[, cell_type_tmp] <- as.numeric(renal.markers$gene %in% marker_genes_tmp)
}
for (biomarker_type_tmp in unique(gene2cellType_tab$Biomarker_Type[gene2cellType_tab$Cell_Type_Abbr == "Other"])) {
  marker_genes_tmp <- gene2cellType_tab$Gene[gene2cellType_tab$Biomarker_Type == biomarker_type_tmp]
  renal.markers[, biomarker_type_tmp] <- as.numeric(renal.markers$gene %in% marker_genes_tmp)
}

write.table(renal.markers, file = paste0(dir_out,  "NonImmune.DEGs.Pos.CellTypeMarkerAnnotated",  ".", run_id, ".txt"), quote = F, sep = "\t", row.names = F)

marker_col_names <- colnames(renal.markers)[!(colnames(renal.markers) %in% c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2"))]
renal.markers.filtered <- renal.markers[rowSums(renal.markers[,marker_col_names]) > 0,]
marker_col_names2print <- marker_col_names[colSums(renal.markers[,marker_col_names]) > 0]
renal.markers.filtered <- renal.markers.filtered[, c("gene", "cluster", marker_col_names2print, "p_val_adj", "avg_logFC", "p_val",  "pct.1", "pct.2")]
write.table(renal.markers.filtered, file = paste0(dir_out,  "NonImmune.DEGs.Pos.CellTypeMarkerOnly", ".", run_id, ".txt"), quote = F, sep = "\t", row.names = F)

###########################################
######## Write # Cells per Cluster across Samples
###########################################
sn_cell_num_tab <- data.frame(renal.int.nonimmune.obj@meta.data %>%
                                select(orig.ident, seurat_clusters) %>%
                                table())
sn_cell_num_tab %>%
  head()
colnames(sn_cell_num_tab) <- c("SampID", "seurat_clusters", "Barcode_Num")
write.table(x = sn_cell_num_tab, file = paste0(dir_out, "Barcodes_per_Cluster_Sample", ".", run_id ,".tsv"), quote = F, sep = "\t", row.names = F)


###########################################
######## Write # Cells per Cluster per Sample
###########################################
num_sn_by_sample <- renal.int.object@meta.data %>% 
  group_by(orig.ident) %>%
  summarize(n())

sn_cell_num_tab <- data.frame(renal.int.nonimmune.obj@meta.data %>%
                                select(orig.ident, seurat_clusters) %>%
                                table())
sn_cell_num_tab %>%
  head()
colnames(sn_cell_num_tab) <- c("SampID", "seurat_clusters", "Barcode_Num")

sn_cell_num_tab <- merge(sn_cell_num_tab, num_sn_by_sample, by.x = c("SampID"), by.y = c("orig.ident"), all.x = T)
colnames(sn_cell_num_tab) <- c("SampID", "seurat_clusters", "Barcode_Num", "Barcode_Sum")

sn_cell_num_tab <- sn_cell_num_tab %>%
  mutate(Barcode_Perc = Barcode_Num/Barcode_Sum)
write.table(x = sn_cell_num_tab, file = paste0(dir_out, "Barcodes_per_Cluster_Sample", ".", run_id ,".tsv"), quote = F, sep = "\t", row.names = F)


