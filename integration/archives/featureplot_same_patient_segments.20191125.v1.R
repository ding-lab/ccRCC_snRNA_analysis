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
case_ids <- c("C3N-01200", "C3N-00733")
# case_ids <- c("C3N-00733")

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191024.v1.tsv")

# plot dotplot by sample --------------------------------------------------
for (case_id_tmp in case_ids) {
  dir_out_tmp <- paste0(dir_out, case_id_tmp, "/")
  dir.create(dir_out_tmp)
  
  if (case_id_tmp == "C3N-00733") {
    seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_segments/20191125.v1/", case_id_tmp, ".Tummor_Segments.Integrated.20191125.v1.RDS")
  }
  if (case_id_tmp == "C3N-01200") {
    seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_segments/20191122.v1/", case_id_tmp, ".Tummor_Segments.Integrated.20191122.v1.RDS")
  }
  
  ## input seurat object
  seurat_obj <- readRDS(file = seurat_obj_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  p <- DimPlot(seurat_obj, reduction = "umap",  label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  
  for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr)) {
    dir_out_cell_type_tmp <- paste0(dir_out_tmp, cell_type_tmp, "/")
    dir.create(dir_out_cell_type_tmp)
    
    genes2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
    genes2plot <- intersect(rownames(seurat_obj@assays$RNA@counts), genes2plot)
    
    for (gene_tmp in genes2plot) {
      ## plot feature plot
      DefaultAssay(seurat_obj) <- "RNA"
      p <- FeaturePlot(object = seurat_obj, features = gene_tmp, 
                       cols = c("grey", "red"), reduction = "umap", label = F)
      p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
      
      p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                              axis.text.y=element_blank(),axis.ticks=element_blank(),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank())
      file2write <- paste(dir_out_cell_type_tmp, case_id_tmp, ".Tumor_Segments.Featureplot.", gene_tmp, "." , run_id, ".png", sep="")
      png(file = file2write, width = 900, height = 800, res = 150)
      print(p)
      dev.off()
    }
    
  }
}


