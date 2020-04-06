# Yige Wu @WashU Oct 2019
## for plotting the marker genes for integrated object, showing cell of origin

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
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191022.v5.tsv")

# input non-immune integrated object --------------------------------------
object2plot <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal_Integrated.20191021.v1.RDS")

# input DEG ---------------------------------------------------------------
deg_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal.DEGs.Pos.txt", data.table = F)

# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191022.v1.tsv", data.table = F)
nonimmune_clusters <- cluster2celltype_tab %>%
  filter(Is_Immune == "No") %>%
  select(Cluster)
nonimmune_clusters <- nonimmune_clusters$Cluster


# Featureplot gene split by sample ----------------------------------------
dir_markers_exp <- paste0(dir_out, "Marker_Expression_by_gene/")
dir.create(dir_markers_exp)
DefaultAssay(object2plot) <- "RNA"

markers2plot <- c("MET")
for (marker_tmp in markers2plot) {
  biomarker_type_tmp <- paste0(gene2cellType_tab$Cell_Type_Abbr[gene2cellType_tab$Gene == marker_tmp], "_", gene2cellType_tab$Biomarker_Type[gene2cellType_tab$Gene == marker_tmp])
  dir_cell_type_markers_exp <- paste0(dir_markers_exp, biomarker_type_tmp, "/")
  dir.create(dir_cell_type_markers_exp)
  
  file2write <- paste0(dir_cell_type_markers_exp, marker_tmp, "_Marker_Expression.", run_id, ".png")
  png(file2write, 
      width = 4000,
      height = 700, res = 150)
  p <- FeaturePlot(object = object2plot, features = marker_tmp, 
                   cols = c("grey", "red"), reduction = "umap", label = T, split.by = "orig.ident", ncol = 4)
  print(p)
  dev.off()
  
}



# Dimplot for markers showing cell of origin ------------------------------
dir_marker_exp <- paste0(dir_out, "Cell_Type_Marker_Expression/")
dir.create(dir_marker_exp)

for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr[gene2cellType_tab$Cell_Type_Abbr != "Other" & gene2cellType_tab$Is_Immune_Cell_Type_Marker == "Non_Immune"])) {
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




# Dimplot for markers by gene ------------------------------
dir_markers_exp <- paste0(dir_out, "Marker_Expression_by_gene/")
dir.create(dir_markers_exp)
DefaultAssay(object2plot) <- "RNA"
for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr[gene2cellType_tab$Cell_Type_Abbr != "Other" & gene2cellType_tab$Is_Immune_Cell_Type_Marker == "Non_Immune"])) {
  markers2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  markers2plot <- markers2plot[markers2plot %in% rownames(object2plot@assays$RNA@counts)]
  
  if (length(markers2plot) > 0) {
    dir_cell_type_markers_exp <- paste0(dir_markers_exp, cell_type_tmp, "/")
    dir.create(dir_cell_type_markers_exp)
    
    for (marker_tmp in markers2plot) {
      file2write <- paste0(dir_cell_type_markers_exp, marker_tmp, "_Marker_Expression.", run_id, ".png")
      png(file2write, 
          width = 1000,
          height = 900, res = 150)
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
      file2write <- paste0(dir_cell_type_markers_exp, marker_tmp, "_Marker_Expression.", run_id, ".png")
      png(file2write, 
          width = 1000,
          height = 900, res = 150)
      p <- FeaturePlot(object = object2plot, features = marker_tmp, 
                       cols = c("grey", "red"), reduction = "umap", label = T)
      print(p)
      dev.off()
    }
  }
}


# Dotplot for all druggable targets -------------------------------------------
DefaultAssay(object2plot) <- "RNA"

genes2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Biomarker_Type %in% c("RCC_targeted_drug", "Immunotherapy")]
genes2plot

p <- DotPlot(object = object2plot, features = genes2plot, col.min = 0)
p$data$gene_biomarker_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Biomarker_Type)
p$data$cluster_cell_type <- plyr::mapvalues(p$data$id, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
p$data$is_immune <- ifelse(p$data$id %in% nonimmune_clusters, "Non_Immune", "Immune")
p$data$is_immune <- factor(p$data$is_immune, levels = c("Non_Immune", "Immune"))

p <- p  + RotatedAxis()
p <- p + facet_grid(is_immune + cluster_cell_type~gene_biomarker_type, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 9),
               strip.placement = "outside")

file2write <- paste0(dir_out, "Dotplot_Immunotherapy_Marker_Exp.", run_id, ".png")
png(file = file2write, width = 2000, height = 1000, res = 150)
print(p)
dev.off()



# plot dimention plot to show the cell types ------------------------------
object2plot@meta.data$cluster_cell_type <- plyr::mapvalues(object2plot@meta.data$seurat_clusters, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
p <- DimPlot(object2plot, reduction = "umap", group.by = "cluster_cell_type", label = T, label.size	= 5, repel = T,)

file2write <- paste(dir_out, "Dimplot_by_Clusters_with_CellType.", run_id, ".png", sep="")
png(file = file2write, width = 2200, height = 1200, res = 150)
print(p)
dev.off()



