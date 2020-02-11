# Yige Wu @WashU Oct 2019
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
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191022.v4.tsv")


# input non-immune integrated object --------------------------------------
object2plot <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/process_nonimmune_cells_on_cluster/20191022.v1/NonImmune_Integrated.20191022.v1.RDS")

# input DEG ---------------------------------------------------------------
deg_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/process_nonimmune_cells_on_cluster/20191022.v1/NonImmune.DEGs.Pos.txt", data.table = F)

# plot DEG in dot plot before cell type assignment-----------------------------------------------------------
DefaultAssay(object2plot) <- "RNA"

genes2plot <- deg_tab$gene %>% unique()
genes2plot <- genes2plot[genes2plot %in% gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]]
genes2plot <- genes2plot[genes2plot %in% gene2cellType_tab$Gene[gene2cellType_tab$Is_Immune_Cell_Type_Marker != "Immune"]]
genes2plot <- c("CA9", "HIF1A", "EPAS1", "VHL", "VIM1", "MKI67", "CASP3", genes2plot)
genes2plot

p <- DotPlot(object = object2plot, features = genes2plot, col.min = 0, dot.min = 0.00000001)
p$data$cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
p <- p  + RotatedAxis()
p <- p + facet_grid(~cell_type, scales = "free_x", space = "free", shrink = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.placement = "outside")
file2write <- paste0(dir_out, "Dotplot_NonImmune_DEG_Marker_Exp_without_CellType.", run_id, ".pdf")
pdf(file2write, width = 15, height = 7)
print(p)
dev.off()

file2write <- paste0(dir_out, "Dotplot_NonImmune_DEG_Marker_Exp_without_CellType.", run_id, ".png")
png(file = file2write, width = 2000, height = 1000, res = 150)
print(p)
dev.off()

# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - NonImmuneCluster2CellType.20191022.v1.tsv", data.table = F)
object2plot@meta.data$cluster_cell_type <- plyr::mapvalues(object2plot@meta.data$seurat_clusters, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)

# plot dot plot with biomarkers other than cell type markers-----------------------------------------------------------
genes2plot <- deg_tab$gene %>% unique()
genes2plot <- genes2plot[genes2plot %in% gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]]
genes2plot <- genes2plot[genes2plot %in% gene2cellType_tab$Gene[gene2cellType_tab$Is_Immune_Cell_Type_Marker != "Immune"]]
genes2plot <- c("CA9", "HIF1A", "EPAS1", "VHL", "VIM1", "MKI67", "CASP3", genes2plot)
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

file2write <- paste0(dir_out, "Dotplot_NonImmune_Marker_Exp.", run_id, ".png")
png(file = file2write, width = 2000, height = 1000, res = 150)
print(p)
dev.off()


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

dir_markers_exp <- paste0(dir_out, "Marker_Expression_by_gene/")
dir.create(dir_markers_exp)


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


