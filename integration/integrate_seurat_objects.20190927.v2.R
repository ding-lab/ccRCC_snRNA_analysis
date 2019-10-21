# Yige Wu @WashU Sep 2019
## for integrating two snRNA datasets for sample CPT0086820004 and CPT0075130004 (from cellranger output with premrna reference)

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# set parameters ----------------------------------------------------------
version_tmp <- 1
sample_ids <- c("CPT0075130004_notFACS", "CPT0086820004_notFACS", "CPT0075140002", "CPT0001260013", "CPT0086350004")
seurat_obj_path_tab <- data.frame(sample_id = sample_ids,
                                  seurat_object_path=paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", sample_ids, ""))
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/scRNA_auto/summary/RCC_snRNA_Downstream_Processing - Seurat_Summary.09272019.v5.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object

###########################################
######## Dataset preprocessing
###########################################
renal.list <- list()
# Input seurat objects -----------------------------------------------------
for (i in 1:nrow(seurat_summary2process)) {
  sample_id_tmp <- seurat_summary$Aliquot[i]
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
reference.list <- renal.list
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


saveRDS(object = renal.integrated, file = paste0(makeOutDir(), run_id, "/renal_integrated", run_id, ".RDS"))
renal.integrated <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/scRNA/intergration/integrate_seurat_objects/20190927.v1/renal_integrated.20190927.v1.RDS")

###########################################
######## Write row read count table
###########################################
raw_exp_df <- renal.integrated@assays$RNA@counts
raw_exp_df <- as.data.frame(raw_exp_df, keep.rownames = T)

fwrite(x = raw_exp_df, file = paste0(dir_out, "raw_counts", ".", run_id ,".tsv"), quote = F, sep = "\t", row.names = T)


###########################################
######## Plotting dimension
###########################################

## make sure the grouping variable is in the meta data
renal.integrated@meta.data %>%
  head()

file2write <- paste(dir_out, "DimPlot_raw_cell_clusters_in_sample_", run_id, ".pdf", sep="")
pdf(file = file2write, width = 8, height = 6)
p2 <- DimPlot(renal.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)
print(p2)
dev.off()

file2write <- paste(dir_out, "DimPlot_by_cluster_and_sample_in_sample_", run_id, ".pdf", sep="")
pdf(file = file2write, width = 15, height = 6)
p1 <- DimPlot(renal.integrated, reduction = "umap", group.by = "orig.ident", label = T)
p2 <- DimPlot(renal.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)
plot_grid(p1, p2)
dev.off()

file2write <- paste(dir_out, "DimPlot_raw_cell_clusters_in_sample_splitted.", run_id, ".pdf", sep="")
pdf(file = file2write, width = 20, height = 6)
DimPlot(renal.integrated, reduction = "umap", split.by = "orig.ident", label = T)
dev.off()

###########################################
######## Differential expression Analysis
###########################################

# find DEG ----------------------------------------------------------------
DefaultAssay(renal.integrated) <- "RNA"

renal.markers <- FindAllMarkers(object = renal.integrated, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
renal.markers %>%
  colnames()
renal.markers <- renal.markers[, c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
write.table(renal.markers, file = paste0(dir_out, "Renal.DEGs.Pos.txt"), quote = F, sep = "\t", row.names = F)

# filter DEG by manual curated markers ------------------------------------
# gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/dependencies/merge_proteinatlas_w_manual_markers/RCC_Marker_Tab_w.HumanProteinAtlast.20190913.v1.tsv")
# gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene_to_Cell_Type_Table.20190916.v2.tsv")
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - RCC_Marker_Tab_w.HumanProteinAtlas.20190925.v1.tsv")

for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr)) {
  marker_genes_tmp <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  renal.markers[, cell_type_tmp] <- as.numeric(renal.markers$gene %in% marker_genes_tmp)
}
write.table(renal.markers, file = paste0(dir_out,  "Renal.DEGs.Pos.CellTypeMarkerAnnotated",  ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp, ".txt"), quote = F, sep = "\t", row.names = F)

renal.markers.filtered <- renal.markers[rowSums(renal.markers[,unique(gene2cellType_tab$Cell_Type_Abbr)]) > 0,]
cell_types2print <- unique(gene2cellType_tab$Cell_Type_Abbr)[colSums(renal.markers[,unique(gene2cellType_tab$Cell_Type_Abbr)]) > 0]
renal.markers.filtered <- renal.markers.filtered[, c("gene", "cluster", cell_types2print, "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
write.table(renal.markers.filtered, file = paste0(dir_out,  "Renal.DEGs.Pos.CellTypeMarkerOnly", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp, ".txt"), quote = F, sep = "\t", row.names = F)

###########################################
######## Plotting marker expression
###########################################
object2plot <- renal.integrated
DefaultAssay(object2plot) <- "RNA"

dir_marker_exp <- paste0(dir_out, "Marker_Expression/")
dir.create(dir_marker_exp)

for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr)) {
  markers2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  file2write <- paste0(dir_marker_exp, cell_type_tmp, "_Marker_in_sample_", run_id, "_FeaturePlot.pdf")
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

for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr)) {
  markers2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  markers2plot <- markers2plot[markers2plot %in% rownames(object2plot@assays$RNA@counts)]
  
  if (length(markers2plot) > 0) {
    dir_cell_type_markers_exp <- paste0(dir_markers_exp, cell_type_tmp, "/")
    dir.create(dir_cell_type_markers_exp)
    
    for (marker_tmp in markers2plot) {
      file2write <- paste0(dir_cell_type_markers_exp, marker_tmp, "_in_sample_", run_id, "_FeaturePlot.pdf")
      
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

###########################################
######## Write # Cells per Cluster
###########################################
sn_cell_num_tab <- data.frame(renal.integrated@meta.data %>%
                                select(orig.ident, seurat_clusters) %>%
                                table())
sn_cell_num_tab %>%
  head()
colnames(sn_cell_num_tab) <- c("SampID", "seurat_clusters", "Barcode_Num")
version_tmp <- 1
write.table(x = sn_cell_num_tab, file = paste0(dir_out, "snRNA_bc_num_per_cluster", ".", run_id ,".tsv"), quote = F, sep = "\t", row.names = F)

###########################################
######## Label Cells
###########################################

# # Label Cells -  --------
# n_clusters <- renal.integrated@meta.data$seurat_clusters %>% nlevels
# current.cluster.ids <- 0:(n_clusters-1)
# new.cluster.ids<- c('MacroMono','MET_ccRCC', "TAM", 'CA9_ccRCC',
#                     'CA9_MET_ccRCC', 'HNF1B_ENPP3_ccRCC','TAM','Mesenchymal',
#                     'CD4Tcell', 'TAM', 'TAM','Endothelial_Myofibroblast',
#                     'Endothelial_Podocyte','Proliferating_TAM',"TAM", "TAM")
# new.cluster.ids <- paste0(new.cluster.ids, "_C", current.cluster.ids)
# renal.integrated@meta.data$cell_type <- plyr::mapvalues(as.character(Idents(renal.integrated)), from= current.cluster.ids,to = new.cluster.ids)
# table(renal.integrated@meta.data$cell_type)
# 
# 
# for (sample_id_tmp in sample_ids) {
#   tab_anno_10xmapping <- renal.integrated@meta.data
#   tab_anno_10xmapping$barcode_integrated <- rownames(tab_anno_10xmapping)
#   tab_anno_10xmapping <- tab_anno_10xmapping %>%
#     mutate(barcode = str_split_fixed(string = barcode_integrated, pattern = "_", n = 2)[,1]) %>%
#     filter(orig.ident == sample_id_tmp) %>%
#     select(barcode, cell_type)
#   
#   write.table(x = tab_anno_10xmapping, file = paste0(makeOutDir(), "tab_anno_10xmapping_", sample_id_tmp, ".tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
# }
# 
# 
# Idents(renal.integrated) <- renal.integrated@meta.data$cell_type
# saveRDS(object = renal.integrated, file = paste0(makeOutDir(), run_id, "_renal_integrated_wCellType.20190923.v1.RDS"))
# 
# 
# 
