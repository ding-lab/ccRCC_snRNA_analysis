# Yige Wu @WashU Jan 2020
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
case_ids <- c("C3N-01200")

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20200131.v1.tsv")

# plot dotplot by sample --------------------------------------------------
for (case_id_tmp in case_ids) {
  if (case_id_tmp == "C3N-00733") {
    seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_tumor_normal/20191125.v1/", case_id_tmp, ".Tumor_Normal.Integrated.20191125.v1.RDS")
    deg_tab_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/run_deg_same_patient_tumor_normal/20191125.v1/", case_id_tmp, ".Tumor_Normal.DEGs.Pos.txt")
  }
  if (case_id_tmp == "C3N-01200") {
    seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_tumor_normal/20200117.v1/", case_id_tmp, ".Tummor_Normal.Integrated.20200117.v1.RDS")
    deg_tab_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/run_deg_same_patient_tumor_normal/20200117.v1/", case_id_tmp, ".Tumor_Normal.DEGs.Pos.txt")
  }
  
  ## input seurat object
  seurat_obj <- readRDS(file = seurat_obj_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  ## input DEG
  deg_tab_path
  deg_tab <- fread(input = deg_tab_path, data.table = F)
  
  ## get marker genes
  marker_genes <- FindAllMarkers(object = seurat_obj, test.use = "roc", only.pos = T, return.thresh = 0.5)
  
  ## examine cell type markers
  cluster_tmp <- 11
  cluster_markers <- marker_genes %>%
    filter(cluster == cluster_tmp) %>%
    filter(gene %in% gene2cellType_tab$Gene)
  
  ## plot marker gene expression
  genes2plot <- gene2cellType_tab$Gene
  genes2plot <- intersect(deg_tab$gene, genes2plot)
  
  p <- DotPlot(object = seurat_obj, features = genes2plot, col.min = 0)
  p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Group)
  p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type1)
  p <- p  + RotatedAxis()
  p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1, scales = "free", space = "free", drop = T)
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 strip.background = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 strip.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
                 axis.text.x = element_text(size = 9),
                 strip.placement = "outside")
  
  file2write <- paste0(dir_out, case_id_tmp,".Tumor_Normal.Dotplot.CellTypeMarkersInDEGs.", run_id, ".png")
  png(file = file2write, width = 2000, height = 1800, res = 150)
  print(p)
  dev.off()
  
  ## plot all marker gene expression
  genes2plot <- gene2cellType_tab$Gene
  # genes2plot <- intersect(rownames(seurat_obj@assays$RNA@counts), genes2plot)
  
  p <- DotPlot(object = seurat_obj, features = genes2plot, col.min = 0)
  p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Group)
  p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type1)
  p <- p  + RotatedAxis()
  p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1, scales = "free", space = "free", drop = T)
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 strip.background = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 strip.text.x = element_text(angle = 90, vjust = 0.5),
                 axis.text.x = element_text(size = 9),
                 strip.placement = "outside")
  
  file2write <- paste0(dir_out,  case_id_tmp,".Tumor_Normal.Dotplot.CellTypeMarkers.", run_id, ".png")
  png(file = file2write, width = 5500, height = 1800, res = 150)
  print(p)
  dev.off()
  
}


