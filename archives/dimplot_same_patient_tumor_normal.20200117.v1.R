# Yige Wu @WashU Nov 2019
## for plotting the marker genes for integrated object per case

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

# set case ids to be processed --------------------------------------------
case_ids <- c("C3N-00733")

# input cluster2celltype table --------------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Tumor_Segments.AllCluster2Cell_Type.20191125.v1.tsv", data.table = F)

# plot dotplot by sample --------------------------------------------------
case_id_tmp <- "C3N-00733"

for (case_id_tmp in case_ids) {
  if (case_id_tmp == "C3N-00733") {
    seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_segments/20191125.v1/", case_id_tmp, ".Tummor_Segments.Integrated.20191125.v1.RDS")
  }
  if (case_id_tmp == "C3N-01200") {
    seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_segments/20191122.v1/", case_id_tmp, ".Tummor_Segments.Integrated.20191122.v1.RDS")
  }
  
  ## input seurat object
  seurat_obj <- readRDS(file = seurat_obj_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  ## get the cluster2celltype info for current aliquot
  cluster2celltype_tab_tmp <- cluster2celltype_tab %>%
    filter(Case == case_id_tmp)

  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj@meta.data$cluster_cell_type <- plyr::mapvalues(seurat_obj@meta.data$seurat_clusters, from = cluster2celltype_tab_tmp$Cluster, to = cluster2celltype_tab_tmp$Enriched_Cell_Type_Abbr)
  seurat_obj@meta.data$cluster_cell_type_uniq <- paste0(seurat_obj@meta.data$cluster_cell_type, "_C", seurat_obj@meta.data$seurat_clusters)
  p <- DimPlot(seurat_obj, reduction = "umap",  group.by = "cluster_cell_type_uniq", label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  
  
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "cluster_cell_type_uniq", label = F, label.size	= 5)
  p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = cluster_cell_type_uniq), box.padding = 0.5)
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          legend.position = "none")
  file2write <- paste(dir_out, case_id_tmp, ".Tumor_Segments.Dimplot.ByCellType.", run_id, ".png", sep="")
  png(file = file2write, width = 1200, height = 1200, res = 150)
  print(p)
  dev.off()
  
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident", label = F, label.size	= 5)
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          legend.position = "top")
  file2write <- paste(dir_out, case_id_tmp, ".Tumor_Segments.Dimplot.BySample.", run_id, ".png", sep="")
  png(file = file2write, width = 1200, height = 1300, res = 150)
  print(p)
  dev.off()
}


