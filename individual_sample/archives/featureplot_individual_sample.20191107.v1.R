# Yige Wu @WashU Nov 2019
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
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191024.v1.tsv")

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
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20191106.v1.tsv", data.table = F)

# input the HIF downstream to plot ----------------------------------------
HIF_targets2plot <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/proteomics/plot_HIF_pahtway_protein_heatmap/ccRCC_snRNA_Downstream_Processing - HIF_Target_Summary.tsv", data.table = F)

# plot featureplot by sample and then by gene--------------------------------------------------
snRNA_aliquot_id_tmp <- "CPT0001220012"
gene_tmp <- "CA9"
gene_tmp <- "HIF1A"
genes2plot <- c("CA9", "VHL", "HIF1A", "EPAS1")
genes2plot <- c("CA9", "VHL", "HIF1A", "EPAS1", "PBRM1", "BAP1", "SETD2")
genes2plot <- HIF_targets2plot$Gene_Symbol[HIF_targets2plot$Tumor == 1 & is.na(HIF_targets2plot$Stromal) & is.na(HIF_targets2plot$Immune)]
genes2plot <- c("ARNT", "ARNT2", genes2plot)

for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  ## get the cluster2celltype info for current aliquot
  cluster2celltype_tab_tmp <- cluster2celltype_tab %>%
    filter(Aliquot == snRNA_aliquot_id_tmp)
  
  seurat_obj@meta.data$cluster_cell_type <- plyr::mapvalues(seurat_obj@meta.data$seurat_clusters, from = cluster2celltype_tab_tmp$Cluster, to = cluster2celltype_tab_tmp$Enriched_Cell_Type_Abbr)
  
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "cluster_cell_type", label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  label_data <- label_data %>%
    filter(cluster_cell_type %in% cluster2celltype_tab_tmp$Enriched_Cell_Type_Abbr[cluster2celltype_tab_tmp$Is_Malignant == "Yes"])
  
  for (gene_tmp in genes2plot) {
    ## plot feature plot
    DefaultAssay(seurat_obj) <- "RNA"
    p <- FeaturePlot(object = seurat_obj, features = gene_tmp, 
                     cols = c("grey", "red"), reduction = "umap", label = F, min.cutoff = "q10", max.cutoff = "q90")
    p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = cluster_cell_type))
    
    p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank())
    file2write <- paste(dir_out, snRNA_aliquot_id_tmp, ".Featureplot.", gene_tmp, ".q10q90." , run_id, ".png", sep="")
    png(file = file2write, width = 900, height = 800, res = 150)
    print(p)
    dev.off()
    
    ## plot feature plot
    DefaultAssay(seurat_obj) <- "RNA"
    p <- FeaturePlot(object = seurat_obj, features = gene_tmp, 
                     cols = c("grey", "red"), reduction = "umap", label = F)
    p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = cluster_cell_type))
    
    p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank())
    file2write <- paste(dir_out, snRNA_aliquot_id_tmp, ".Featureplot.", gene_tmp, "." , run_id, ".png", sep="")
    png(file = file2write, width = 900, height = 800, res = 150)
    print(p)
    dev.off()
  }
}



# plot featureplot by sample and then by VHLX(HIF targets)--------------------------------------------------
snRNA_aliquot_id_tmp <- "CPT0001220012"
gene_tmp <- "CA9"
genes2plot <- c("CA9", "VHL", "HIF1A", "EPAS1", "PBRM1", "BAP1", "SETD2")

for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  ## get the cluster2celltype info for current aliquot
  cluster2celltype_tab_tmp <- cluster2celltype_tab %>%
    filter(Aliquot == snRNA_aliquot_id_tmp)
  
  seurat_obj@meta.data$cluster_cell_type <- plyr::mapvalues(seurat_obj@meta.data$seurat_clusters, from = cluster2celltype_tab_tmp$Cluster, to = cluster2celltype_tab_tmp$Enriched_Cell_Type_Abbr)
  
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "cluster_cell_type", label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  label_data <- label_data %>%
    filter(cluster_cell_type %in% cluster2celltype_tab_tmp$Enriched_Cell_Type_Abbr[cluster2celltype_tab_tmp$Is_Malignant == "Yes"])
  
  for (gene_tmp in genes2plot) {
    
    ## plot feature plot
    DefaultAssay(seurat_obj) <- "RNA"
    p <- FeaturePlot(object = seurat_obj, features = c("VHL", gene_tmp), blend = T, reduction = "umap", label = F)
    # p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = cluster_cell_type))
    p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank())
    file2write <- paste(dir_out, snRNA_aliquot_id_tmp, ".Featureplot.", gene_tmp, "." , run_id, ".png", sep="")
    png(file = file2write, width = 900, height = 800, res = 150)
    print(p)
    dev.off()
    
    p <- FeatureScatter(seurat_obj, feature1 = "VHL", feature2 = gene_tmp)
    p
  }
}

