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

# plot dotplot by sample --------------------------------------------------
snRNA_aliquot_id_tmp <- "CPT0001220012"
marker2plot <- "CA9"
marker2plot <- "HIF1A"
markers2plot <- c("CA9", "VHL", "HIF1A", "EPAS1")
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
    for (marker2plot in markers2plot) {
      
    ## input seurat object
    seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
    seurat_obj <- readRDS(file = seurat_obj_path)
    
    ## get the cluster2celltype info for current aliquot
    cluster2celltype_tab_tmp <- cluster2celltype_tab %>%
      filter(Aliquot == snRNA_aliquot_id_tmp)
    
    DefaultAssay(seurat_obj) <- "RNA"
    p <- FeaturePlot(object = seurat_obj, features = marker2plot, 
                     cols = c("grey", "red"), reduction = "umap", label = F, min.cutoff = "q10", max.cutoff = "q90")
    p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank())
    file2write <- paste(dir_out, snRNA_aliquot_id_tmp, ".Featureplot.", marker2plot, ".", run_id, ".png", sep="")
    png(file = file2write, width = 900, height = 800, res = 150)
    print(p)
    dev.off()
  }
}
