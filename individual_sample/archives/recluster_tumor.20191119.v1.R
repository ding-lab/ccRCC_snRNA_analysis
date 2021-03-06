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

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191021.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds")) %>%
  mutate(Paht_deg_table = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                 "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                 "/", Aliquot, FACS, ".DEGs.Pos.txt"))
seurat_summary2process$Path_seurat_object

# input cluster2celltype table --------------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20191119.v2.tsv", data.table = F)


# set aliquot ids to be processed -----------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# set aliquot id ----------------------------------------------------------
snRNA_aliquot_id_tmp <- "CPT0019130004"

for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  # set aliquot specific output directory -----------------------------------
  dir_out_tmp <- paste0(dir_out, snRNA_aliquot_id_tmp, "/")
  dir.create(dir_out_tmp)
  
  # Input the seurat object -------------------------------
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj_path
  orig.object <- readRDS(file = seurat_obj_path)
  
  
  # set clusters to be processed --------------------------------------------
  clusters2process <- cluster2celltype_tab %>%
    filter(Aliquot == snRNA_aliquot_id_tmp) %>%
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
  saveRDS(object = new.object, file = paste0(dir_out_tmp, snRNA_aliquot_id_tmp, ".Malignant_Reclustered.", run_id, ".RDS"))
  
  # get DEG ----------------------------------------------------------------
  DefaultAssay(new.object) <- "RNA"
  
  deg_tab <- FindAllMarkers(object = new.object, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  deg_tab %>%
    head()
  deg_tab <- deg_tab[, c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
  write.table(deg_tab, file = paste0(dir_out_tmp, snRNA_aliquot_id_tmp, ".Malignant_Reclustered", ".DEGs.Pos.txt"), quote = F, sep = "\t", row.names = F)
  
  # plot dimention reduction ---------------------------------------------
  file2write <- paste(dir_out_tmp, snRNA_aliquot_id_tmp, ".Malignant_Reclustered.Dimplot.", run_id, ".png", sep="")
  png(file = file2write, width = 850, height = 800, res = 150)
  p <- DimPlot(new.object, reduction = "umap", group.by = "seurat_clusters", label = T)
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  
  print(p)
  dev.off()
  
  
}
 
