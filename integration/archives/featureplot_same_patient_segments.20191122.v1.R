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
case_ids <- c("C3N-01200", "C3L-01313")

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191024.v1.tsv")

# plot dotplot by sample --------------------------------------------------
case_id_tmp <- "C3L-01313"
for (case_id_tmp in case_ids) {
  dir_out_tmp <- paste0(dir_out, case_id_tmp, "/")
  dir.create(dir_out_tmp)
  
  ## input seurat object
  seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_segments/20191122.v1/", case_id_tmp, ".Tummor_Segments.Integrated.20191122.v1.RDS")
  seurat_obj <- readRDS(file = seurat_obj_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  ## input DEG
  deg_tab_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/run_deg_same_patient_segments/20191122.v1/", case_id_tmp, ".Tumor_Segments.DEGs.Pos.txt")
  deg_tab_path
  deg_tab <- fread(input = deg_tab_path, data.table = F)
  
  ## plot marker gene expression
  genes2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]
  genes2plot <- intersect(deg_tab$gene, genes2plot)
  genes2plot <- c("CA9", genes2plot)
  
  p <- DimPlot(seurat_obj, reduction = "umap",  label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  
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
    file2write <- paste(dir_out_tmp, case_id_tmp, ".Tumor_Segments.Featureplot.", gene_tmp, "." , run_id, ".png", sep="")
    png(file = file2write, width = 900, height = 800, res = 150)
    print(p)
    dev.off()
  }
}


