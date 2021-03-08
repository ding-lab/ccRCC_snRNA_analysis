# Yige Wu @WashU Sep 2019
## for isolating the non-immune cell clusters and re-do clustering

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# set case ids to be processed --------------------------------------------
case_ids <- c("C3N-00733")

# input cluster2celltype table --------------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Tumor_Segments.AllCluster2Cell_Type.20191125.v1.tsv", data.table = F)

# set aliquot id ----------------------------------------------------------
for (case_id_tmp in case_ids) {
  dir_out_tmp <- paste0(dir_out, case_id_tmp, "/")
  dir.create(dir_out_tmp)
  
  if (case_id_tmp == "C3N-00733") {
    seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_segments/20191125.v1/", case_id_tmp, ".Tummor_Segments.Integrated.20191125.v1.RDS")
  }
  if (case_id_tmp == "C3N-01200") {
    seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_segments/20191122.v1/", case_id_tmp, ".Tummor_Segments.Integrated.20191122.v1.RDS")
  }
  
  # Input the seurat object -------------------------------
  orig.object <- readRDS(file = seurat_obj_path)
  
  # set clusters to be processed --------------------------------------------
  clusters2process <- cluster2celltype_tab %>%
    filter(Case == case_id_tmp) %>%
    filter(Is_Malignant == "Yes") %>%
    select(Cluster)
  
  clusters2process <- clusters2process$Cluster
  clusters2process
  
  # subset object by non-immune clusters ------------------------------------
  new.object <- subset(orig.object, idents = clusters2process)
  
  # Run the standard workflow for visualization and clustering ------------
  new.object <- FindVariableFeatures(object = new.object, selection.method = "vst", nfeatures = 2000)
  new.object <- ScaleData(new.object, verbose = F) 
  new.object <- RunPCA(new.object, npcs = 30, verbose = FALSE)
  new.object <- RunUMAP(new.object, reduction = "pca", dims = 1:20)
  new.object <- FindNeighbors(new.object, reduction = "pca", dims = 1:20)
  new.object <- FindClusters(new.object, resolution = 0.5)
  saveRDS(object = new.object, file = paste0(dir_out_tmp, case_id_tmp, ".Tumor_Segments.Malignant_Reclustered.", run_id, ".RDS"))
  
  # find DEG ----------------------------------------------------------------
  DefaultAssay(new.object) <- "RNA"
  deg_tab <- FindAllMarkers(object = new.object, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  deg_tab %>%
    head()
  deg_tab <- deg_tab[, c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
  write.table(deg_tab, file = paste0(dir_out_tmp, case_id_tmp, ".Tumor_Segments.Malignant_Reclustered", ".DEGs.Pos.txt"), quote = F, sep = "\t", row.names = F)
  
  # plot dimention reduction ---------------------------------------------
  file2write <- paste(dir_out_tmp, case_id_tmp, ".Tumor_Segments.Malignant_Reclustered.Dimplot.", run_id, ".png", sep="")
  png(file = file2write, width = 800, height = 800, res = 150)
  p <- DimPlot(new.object, reduction = "umap", group.by = "seurat_clusters", label = T)
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank(), legend.position = "none")
  
  print(p)
  dev.off()
  
  file2write <- paste(dir_out_tmp, case_id_tmp, ".Tumor_Segments.Malignant_Reclustered.Dimplot.BySample.", run_id, ".png", sep="")
  png(file = file2write, width = 800, height = 850, res = 150)
  p <- DimPlot(new.object, reduction = "umap", group.by = "orig.ident", label = F)
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank(), legend.position = "top")
  
  print(p)
  dev.off()
}
 
