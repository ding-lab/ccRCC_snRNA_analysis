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
  
  p <- DotPlot(object = seurat_obj, features = genes2plot, col.min = 0)
  p$data$gene_cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
  p <- p  + RotatedAxis()
  p <- p + facet_grid(.~gene_cell_type, scales = "free", space = "free", drop = T)
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 strip.background = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 strip.text.x = element_text(angle = 90, vjust = 0.5),
                 axis.text.x = element_text(size = 9),
                 strip.placement = "outside")
  
  file2write <- paste0(dir_out, case_id_tmp,".Tumor_Segments.Dotplot.CellTypeMarkersInDEGs.", run_id, ".png")
  png(file = file2write, width = 2000, height = 1000, res = 150)
  print(p)
  dev.off()
  
  ## plot all marker gene expression
  genes2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]
  genes2plot <- intersect(rownames(seurat_obj@assays$RNA@counts), genes2plot)
  genes2plot <- c("CA9", genes2plot)
  
  p <- DotPlot(object = seurat_obj, features = genes2plot, col.min = 0)
  p$data$gene_cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
  p <- p  + RotatedAxis()
  p <- p + facet_grid(.~gene_cell_type, scales = "free", space = "free", drop = T)
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 strip.background = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 strip.text.x = element_text(angle = 90, vjust = 0.5),
                 axis.text.x = element_text(size = 9),
                 strip.placement = "outside")
  
  file2write <- paste0(dir_out,  case_id_tmp,".Tumor_Segments.Dotplot.CellTypeMarkers.", run_id, ".png")
  png(file = file2write, width = 6000, height = 1000, res = 150)
  print(p)
  dev.off()
}


