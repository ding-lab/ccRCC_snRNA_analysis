# Yige Wu @WashU Nov 2019
## for plotting the marker genes for individual sample in dotplot

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

# plot dotplot by sample --------------------------------------------------
snRNA_aliquot_id_tmp <- "CPT0001180011"
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  ## input DEG
  deg_tab_path <- seurat_summary2process$Paht_deg_table[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
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
  
  file2write <- paste0(dir_out, snRNA_aliquot_id_tmp,".Individual_Clustered.Dotplot.CellTypeMarkersInDEGs.", run_id, ".png")
  png(file = file2write, width = 2000, height = 1000, res = 150)
  print(p)
  dev.off()
  
  ## plot marker gene expression
  genes2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]
  genes2plot <- intersect(rownames(seurat_obj@assays$RNA@scale.data), genes2plot)
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
  
  file2write <- paste0(dir_out, snRNA_aliquot_id_tmp,".Individual_Clustered.Dotplot.CellTypeMarkers.", run_id, ".png")
  png(file = file2write, width = 3000, height = 1000, res = 150)
  print(p)
  dev.off()
  
}


