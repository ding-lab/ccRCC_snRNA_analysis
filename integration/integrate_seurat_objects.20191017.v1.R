# Yige Wu @WashU Sep 2019
## for integrating two snRNA datasets for sample CPT0086820004 and CPT0075130004 (from cellranger output with premrna reference)

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# set aliquot ids to be processed -----------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191015.v2.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object


# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191004.v1.tsv")


# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191017.v1.tsv", data.table = F)


###########################################
######## Dataset preprocessing
###########################################
renal.list <- list()
# Input seurat objects -----------------------------------------------------
for (i in 1:nrow(seurat_summary2process)) {
  sample_id_tmp <- seurat_summary2process$Aliquot[i]
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[i]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  seurat_obj$orig.ident  <- sample_id_tmp
  DefaultAssay(seurat_obj) <- "RNA" 
  seurat_obj@assays$SCT <- NULL
  seurat_obj@graphs <- list()
  seurat_obj@neighbors <- list()
  seurat_obj@reductions <- list()
  for (col_name_tmp in c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.5", "seurat_clusters")) {
    if (col_name_tmp %in% names(seurat_obj@meta.data)) {
      seurat_obj@meta.data[[col_name_tmp]] <- NULL
    }
  }
  renal.list[[sample_id_tmp]] <- seurat_obj
  
}
length(renal.list)

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

saveRDS(object = renal.integrated, file = paste0(dir_out, "Renal_Integrated.", run_id, ".RDS"))

###########################################
######## Plotting dimension
###########################################

## make sure the grouping variable is in the meta data
renal.integrated@meta.data %>%
  select(orig.ident) %>%
  unique()

# plot dimensions by cluster with all aliquots together ------------------------------
p <- DimPlot(renal.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)

file2write <- paste(dir_out, "Dimplot_by_Clusters.", run_id, ".pdf", sep="")
pdf(file = file2write, width = 7, height = 6)
print(p)
dev.off()

file2write <- paste(dir_out, "Dimplot_by_Clusters.", run_id, ".png", sep="")
png(file = file2write, width = 1000, height = 800, res = 150)
print(p)
dev.off()

# plot dimensions by aliquot with all aliquots together ------------------------------
p <- DimPlot(renal.integrated, reduction = "umap", group.by = "orig.ident", label = F)

file2write <- paste(dir_out, "Dimplot_by_Aliquots.", run_id, ".pdf", sep="")
pdf(file = file2write, width = 8, height = 6)
print(p)
dev.off()

file2write <- paste(dir_out, "Dimplot_by_Aliquots.", run_id, ".png", sep="")
png(file = file2write, width = 1100, height = 800, res = 150)
print(p)
dev.off()

# plot dimensions by aliquot and by cluster ------------------------------
p <- DimPlot(renal.integrated, reduction = "umap", label = T, order = snRNA_aliquot_ids, split.by = "orig.ident", ncol = 4)

file2write <- paste(dir_out, "Dimplot_by_Clusters_and_Aliquots.", run_id, ".pdf", sep="")
pdf(file = file2write, width = 20, height = 10)
print(p)
dev.off()

file2write <- paste(dir_out, "Dimplot_by_Clusters_and_Aliquots.", run_id, ".png", sep="")
png(file = file2write, width = 3000, height = 1600, res = 150)
print(p)
dev.off()

###########################################
######## Differential expression Analysis
###########################################

DefaultAssay(renal.integrated) <- "RNA"

renal.markers <- FindAllMarkers(object = renal.integrated, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
renal.markers %>%
  colnames()
renal.markers <- renal.markers[, c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
write.table(renal.markers, file = paste0(dir_out, "Renal.DEGs.Pos.txt"), quote = F, sep = "\t", row.names = F)


list_DEGs_by_cluster <- list()
for (i in unique(renal.markers$cluster)) {
  df2write <- renal.markers %>%
    filter(cluster == i) %>%
    arrange(desc(avg_logFC))
  list_DEGs_by_cluster[[i]] <- df2write
}
file2write <- paste0(dir_out, "Renal.DEGs.Pos.", run_id, ".xlsx")
write.xlsx(list_DEGs_by_cluster, file = file2write)

# filter DEG by manual curated markers ------------------------------------
for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr[gene2cellType_tab$Cell_Type_Abbr != "Other"])) {
  marker_genes_tmp <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  renal.markers[, cell_type_tmp] <- as.numeric(renal.markers$gene %in% marker_genes_tmp)
}
for (biomarker_type_tmp in unique(gene2cellType_tab$Biomarker_Type[gene2cellType_tab$Cell_Type_Abbr == "Other"])) {
  marker_genes_tmp <- gene2cellType_tab$Gene[gene2cellType_tab$Biomarker_Type == biomarker_type_tmp]
  renal.markers[, biomarker_type_tmp] <- as.numeric(renal.markers$gene %in% marker_genes_tmp)
}
write.table(renal.markers, file = paste0(dir_out,  "Renal.DEGs.Pos.CellTypeMarkerAnnotated",  ".", run_id, ".txt"), quote = F, sep = "\t", row.names = F)

marker_col_names <- colnames(renal.markers)[!(colnames(renal.markers) %in% c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2"))]
renal.markers.filtered <- renal.markers[rowSums(renal.markers[,marker_col_names]) > 0,]
marker_col_names2print <- marker_col_names[colSums(renal.markers[,marker_col_names]) > 0]
renal.markers.filtered <- renal.markers.filtered[, c("gene", "cluster", marker_col_names2print, "p_val_adj", "avg_logFC", "p_val",  "pct.1", "pct.2")]
write.table(renal.markers.filtered, file = paste0(dir_out,  "Renal.DEGs.Pos.CellTypeMarkerOnly", ".", run_id, ".txt"), quote = F, sep = "\t", row.names = F)


# plot dot plot before cell type assignment-----------------------------------------------------------
object2plot <- renal.integrated
DefaultAssay(object2plot) <- "RNA"

genes2plot <- renal.markers$gene %>% unique()
genes2plot <- genes2plot[genes2plot %in% gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]]
genes2plot <- c("CA9", "HIF1A", "EPAS1", "VHL", "VIM1", "MKI67", "CASP3", genes2plot)
genes2plot

p <- DotPlot(object = renal.int.nonimmune.obj, features = genes2plot, col.min = 0, dot.min = 0.00000001)
p$data$cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
p <- p  + RotatedAxis()
p <- p + facet_grid(.~cell_type, scales = "free_x", space = "free", shrink = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.placement = "outside")
file2write <- paste0(dir_out, "Dotplot_AllClusters_Marker_Exp_without_CellType.", run_id, ".pdf")
pdf(file2write, width = 16, height = 6)
print(p)
dev.off()

# plot dot plot after cell type assignment-----------------------------------------------------------
object2plot <- renal.integrated
DefaultAssay(object2plot) <- "RNA"

nonimmune_clusters <- cluster2celltype_tab %>%
  filter(Is_Immune == "No") %>%
  select(Cluster)
nonimmune_clusters <- nonimmune_clusters$Cluster

genes2plot <- renal.markers$gene %>% unique()
genes2plot <- genes2plot[genes2plot %in% gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]]
genes2plot <- c("CA9", genes2plot)
genes2plot

p <- DotPlot(object = object2plot, features = genes2plot, col.min = 0)
p$data$gene_cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
p$data$cluster_cell_type <- plyr::mapvalues(p$data$id, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
p$data$gene_cell_type[p$data$features.plot == "CASP3"] <- "Other"
p$data$is_immune <- ifelse(p$data$id %in% nonimmune_clusters, "Non_Immune", "Immune")
p$data$is_immune <- factor(p$data$is_immune, levels = c("Non_Immune", "Immune"))
p <- p  + RotatedAxis()
p <- p + facet_grid(is_immune + cluster_cell_type~gene_cell_type, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 9),
               strip.placement = "outside")

file2write <- paste0(dir_out, "Dotplot_AllClusters_Marker_Exp.", run_id, ".png")
png(file2write, width = 3000, height = 1000, res = 150)
print(p)
dev.off()

# plot dot plot after cell type assignment just for tumor clusters-----------------------------------------------------------
object2plot <- renal.integrated
object2plot <- subset(object2plot, subset = seurat_clusters %in% cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Malignant == "Yes"])
DefaultAssay(object2plot) <- "RNA"

genes2plot <- renal.markers$gene[renal.markers$cluster %in% cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Malignant == "Yes"]]
genes2plot <- intersect(genes2plot, gene2cellType_tab$Gene[!(gene2cellType_tab$Cell_Type_Abbr %in% c("Monocyte", "NK"))])
genes2plot

p <- DotPlot(object = object2plot, features = genes2plot, col.min = 0)
p$data$gene_cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
p$data$cluster_cell_type <- plyr::mapvalues(p$data$id, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
p$data$gene_cell_type[p$data$features.plot == "CASP3"] <- "Other"
p <- p  + RotatedAxis()
p <- p + facet_grid(cluster_cell_type~gene_cell_type, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 9),
               strip.placement = "outside")

file2write <- paste0(dir_out, "Dotplot_TumprClusters_DEG_Exp.", run_id, ".png")
png(file2write, width = 1500, height = 700, res = 150)
print(p)
dev.off()

# plot dimensions after cell type assignment by cluster with all aliquots together ------------------------------
object2plot <- renal.integrated
object2plot@meta.data$cluster_cell_type <- plyr::mapvalues(object2plot@meta.data$seurat_clusters, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
p <- DimPlot(object2plot, reduction = "umap", group.by = "cluster_cell_type", label = T)

file2write <- paste(dir_out, "Dimplot_by_Clusters_with_CellType.", run_id, ".png", sep="")
png(file = file2write, width = 1200, height = 900, res = 150)
print(p)
dev.off()

# plot dimensions after cell type assignment by aliquot and by cluster ------------------------------
p <- DimPlot(object2plot, reduction = "umap", group.by = "cluster_cell_type", label = T, order = snRNA_aliquot_ids, split.by = "orig.ident", ncol = 4)

file2write <- paste(dir_out, "Dimplot_by_Clusters_and_Aliquots_with_CellType.", run_id, ".pdf", sep="")
pdf(file = file2write, width = 20, height = 10)
print(p)
dev.off()

file2write <- paste(dir_out, "Dimplot_by_Clusters_and_Aliquots_with_CellType.", run_id, ".png", sep="")
png(file = file2write, width = 3000, height = 1600, res = 150)
print(p)
dev.off()

p <- DimPlot(object2plot, reduction = "umap", group.by = "cluster_cell_type", label = T, label.size	= 5, repel = T,
             order = snRNA_aliquot_ids, split.by = "orig.ident", ncol = 4)
file2write <- paste(dir_out, "Dimplot_by_Clusters_and_Aliquots_with_CellType.FontAdjusted", run_id, ".png", sep="")
png(file = file2write, width = 3000, height = 1600, res = 150)
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



# Write # Cells per Cluster -----------------------------------------------
sn_cell_num_tab <- data.frame(renal.integrated@meta.data %>%
                                select(orig.ident, seurat_clusters) %>%
                                table())
sn_cell_num_tab %>%
  head()
colnames(sn_cell_num_tab) <- c("snRNA_Aliquot_ID", "Cluster", "Num_Cluster_Barcode")
sn_cell_num_tab$Cluster_Cell_Type <- mapvalues(x = sn_cell_num_tab$Cluster, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
sn_cell_sum_tab <- sn_cell_num_tab %>%
  group_by(snRNA_Aliquot_ID) %>%
  summarise(Num_Aliquot_Barcode = sum(Num_Cluster_Barcode))
sn_cell_num_tab <- merge(sn_cell_num_tab, sn_cell_sum_tab, by = c("snRNA_Aliquot_ID"), all.x = T)
sn_cell_num_tab <- sn_cell_num_tab %>%
  mutate(Perc_Cluster_Barcode = Num_Cluster_Barcode/Num_Aliquot_Barcode)

write.table(x = sn_cell_num_tab, file = paste0(dir_out, "snRNA_Barcode_Num_Per_Cluster_Long", ".", run_id ,".tsv"), quote = F, sep = "\t", row.names = F)

sn_cell_num_mat <- dcast(data = sn_cell_num_tab, formula = snRNA_Aliquot_ID ~ Cluster_Cell_Type, value.var = "Perc_Cluster_Barcode")
write.table(x = sn_cell_num_mat, file = paste0(dir_out, "snRNA_Barcode_Num_Per_Cluster_Wide", ".", run_id ,".tsv"), quote = F, sep = "\t", row.names = F)

###########################################
######## Write annotation file for 10XMapping and inferCNV for each sample and together
###########################################
for (i in 1:nrow(seurat_summary2process)) {
  sample_id_tmp <- seurat_summary2process$Aliquot[i]
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[i]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  tab_anno <- seurat_obj@meta.data %>%
    rownames_to_column("barcode") %>%
    select(barcode, seurat_clusters)
  
}

table2w <- renal.integrated@meta.data
table2w$barcode <- rownames(table2w)
table2w <- table2w %>%
  select(barcode, seurat_clusters)
write.table(x = table2w, file = paste0(dir_out, "snRNA_bc_to_cluster", ".", run_id ,".tsv"), quote = F, sep = "\t", row.names = F, col.names = F)

