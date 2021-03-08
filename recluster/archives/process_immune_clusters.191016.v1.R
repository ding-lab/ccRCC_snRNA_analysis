# Yige Wu @WashU Sep 2019
## for isolating the immune cell clusters and re-do clustering

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input integrated data ---------------------------------------------------
renal.immune.object <- readRDS("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/integration/integrate_seurat_objects/20191015.v1/Renal_Integrated.20191015.v1.RDS")


# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191016.v3.tsv", data.table = F)

# Subset seurat object ----------------------------------------------------
renal.immune.object <- subset(renal.immune.object, subset = seurat_clusters %in% cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Immune == "Yes"] )
rm(renal.immune.object)
snRNA_aliquot_ids <- unique(renal.immune.object@meta.data$orig.ident)
snRNA_aliquot_ids

# Redo integration by creating the list of objects -----------------------------------------------------
renal.list <- list()
for (sample_id_tmp in snRNA_aliquot_ids) {
  seurat_obj_tmp <- subset(renal.immune.object, subset =  orig.ident == sample_id_tmp)
  
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

# Integrate Data -----------------------------------------------------------
renal.anchors <- FindIntegrationAnchors(object.list = renal.list, dims = 1:20)
renal.int.immune.obj <- IntegrateData(anchorset = renal.anchors, dims = 1:20)

# switch to integrated assay. ---------------------------------------------
DefaultAssay(renal.int.immune.obj) <- "integrated" #only have 3000 features

# Run the standard workflow for visualization and clustering ------------
# The variable features of this assay are automatically set during IntegrateData
renal.int.immune.obj <- ScaleData(renal.int.immune.obj, verbose = F) 
renal.int.immune.obj <- RunPCA(renal.int.immune.obj, npcs = 30, verbose = FALSE)
renal.int.immune.obj <- RunUMAP(renal.int.immune.obj, reduction = "pca", dims = 1:20)
renal.int.immune.obj <- FindNeighbors(renal.int.immune.obj, reduction = "pca", dims = 1:20)
renal.int.immune.obj <- FindClusters(renal.int.immune.obj, resolution = 0.5)
saveRDS(object = renal.int.immune.obj, file = paste0(dir_out, "Immune_Integrated.", run_id, ".RDS"))

# find DEG ----------------------------------------------------------------
DefaultAssay(renal.int.immune.obj) <- "RNA"

renal.markers <- FindAllMarkers(object = renal.int.immune.obj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
renal.markers %>%
  colnames()
renal.markers <- renal.markers[, c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
write.table(renal.markers, file = paste0(dir_out, "Immune.DEGs.Pos.txt"), quote = F, sep = "\t", row.names = F)

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191004.v1.tsv")

# plot dot plot before cell type assignment-----------------------------------------------------------
genes2plot <- renal.markers$gene %>% unique()
genes2plot <- genes2plot[genes2plot %in% gene2cellType_tab$Gene[!(gene2cellType_tab$Cell_Type_Abbr %in% c("Other", "PT", "DCT", "MC", "vSMC-P", "Fib", "LOH", "ALH", "DCT_CNT", "ALH_CNT", "CNT", "EC", "Interstitium"))]]
genes2plot

p <- DotPlot(object = renal.int.immune.obj, features = genes2plot, col.min = 0, dot.min = 0.00000001)
p$data$cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
p <- p  + RotatedAxis()
p <- p + facet_grid(.~cell_type, scales = "free_x", space = "free", shrink = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.placement = "outside")
file2write <- paste0(dir_out, "Dotplot_ImmuneClusters_Marker_Exp_without_CellType.", run_id, ".png")
png(file = file2write, width = 2000, height = 1000, res = 150)
print(p)
dev.off()


# plot marker expression in dimention plot --------------------------------
dir_marker_exp <- paste0(dir_out, "Cell_Type_Marker_Expression/")
dir.create(dir_marker_exp)

for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr[gene2cellType_tab$Cell_Type_Abbr != "Other"])) {
  markers2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  file2write <- paste0(dir_marker_exp, cell_type_tmp, "_Marker_Expression.", run_id, ".pdf")
  pdf(file2write, 
      width = 12,
      height = 16)
  p <- FeaturePlot(object = renal.int.immune.obj, features = markers2plot, 
                   cols = c("grey", "red"), reduction = "umap", label = T)
  print(p)
  dev.off()
}

# input cluster 2 cell type table -----------------------------------------
immunecluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - ImmuneCluster2Cell_Type.20191016.v1.tsv", data.table = F)

# plot dimensions after cell type assignment by aliquot and by cluster ------------------------------
renal.int.immune.obj@meta.data$cluster_cell_type <- plyr::mapvalues(renal.int.immune.obj@meta.data$seurat_clusters, from = immunecluster2celltype_tab$Cluster, to = immunecluster2celltype_tab$Enriched_Cell_Type_Abbr)
snRNA_aliquot_ids_w_immune_subtypes <- c("CPT0001180011", "CPT0019130004", "CPT0010110013", "CPT0001260013", "CPT0086350004")
object2plot <- subset(renal.int.immune.obj, subset = orig.ident %in% snRNA_aliquot_ids_w_immune_subtypes)

p <- DimPlot(object2plot, reduction = "umap", group.by = "cluster_cell_type", 
             label = T, label.size	= 4.5, repel = T,
             split.by = "orig.ident", ncol = 3)
file2write <- paste(dir_out, "Dimplot_ImmuneClusters_by_Clusters_and_Aliquots_with_CellType", run_id, ".png", sep="")
png(file = file2write, width = 3000, height = 1600, res = 150)
print(p)
dev.off()


