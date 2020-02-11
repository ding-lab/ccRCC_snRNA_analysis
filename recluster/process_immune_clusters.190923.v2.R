# Yige Wu @WashU Sep 2019
## for isolating the immune cell clusters and re-do clustering

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# set parameters ----------------------------------------------------------
version_tmp <- 3

# Input the identify-assigned seurat object -------------------------------
renal.immune.object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/scRNA/intergration/integrate_seurat_objects_20190914_v1/CPT0075140002_CPT0075130004_CPT0086820004_renal_integrated_wCellType.20190923.v1.RDS")

# Subset seurat object ----------------------------------------------------
renal.immune.object@meta.data$cell_type %>% head()
idents2subset <- c('TAM_C0', "TAM_C2", 'TAM_C6',
                   'CD4Tcell_C8', 'TAM_C9', 'TAM_C10',
                   'Proliferating_TAM_C13',"TAM_C14", "TAM_C15")

renal.immune.object <- subset(renal.immune.object, idents = idents2subset)

# Re-plot dimention reduction ---------------------------------------------
file2write <- paste(makeOutDir(), "DimPlot_immune_cell_clusters_in_sample_", run_id,  ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp, ".pdf", sep="")
pdf(file = file2write, width = 8, height = 6)
p2 <- DimPlot(renal.immune.object, reduction = "umap", group.by = "seurat_clusters", label = T)
print(p2)
dev.off()

###########################################
######## Differential expression
###########################################

# find DEG ----------------------------------------------------------------
DefaultAssay(renal.immune.object) <- "RNA"

immune.markers <- FindAllMarkers(object = renal.immune.object, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
immune.markers %>%
  colnames()
immune.markers <- immune.markers[, c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
write.table(immune.markers, file = paste0(makeOutDir(), "Immune.DEGs.Pos", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp, ".txt"), quote = F, sep = "\t", row.names = F)

# filter DEG by manual curated markers ------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - RCC_Marker_Tab_w.HumanProteinAtlas.20190923.v3.tsv")
immune.markers <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/scRNA/demo/redo_immune_clusters_190923_v1/Immune.DEGs.Pos.txt", data.table = F)

for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr)) {
  marker_genes_tmp <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  immune.markers[, cell_type_tmp] <- as.numeric(immune.markers$gene %in% marker_genes_tmp)
}
write.table(immune.markers, file = paste0(makeOutDir(), "Immune.DEGs.Pos.CellTypeMarkerAnnotated",  ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp, ".txt"), quote = F, sep = "\t", row.names = F)

immune.markers.filtered <- immune.markers[rowSums(immune.markers[,unique(gene2cellType_tab$Cell_Type_Abbr)]) > 0,]
cell_types2print <- unique(gene2cellType_tab$Cell_Type_Abbr)[colSums(immune.markers[,unique(gene2cellType_tab$Cell_Type_Abbr)]) > 0]
immune.markers.filtered <- immune.markers.filtered[, c("gene", "cluster", cell_types2print, "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
write.table(immune.markers.filtered, file = paste0(makeOutDir(), "Immune.DEGs.Pos.CellTypeMarkerOnly", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp, ".txt"), quote = F, sep = "\t", row.names = F)

###########################################
######## Plotting marker expression
###########################################
object2plot <- renal.immune.object
DefaultAssay(object2plot) <- "RNA"

marker_exp_out_path <- paste0(makeOutDir(), "Marker_Expression", "_", format(Sys.Date(), "%Y%m%d") , "_v", version_tmp, "/")
dir.create(marker_exp_out_path)

for (cell_type_tmp in unique(gene2cellType_tab$Cell_Type_Abbr)) {
  markers2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  file2write <- paste0(marker_exp_out_path, cell_type_tmp, "_Marker_in_sample_", run_id, "_FeaturePlot.pdf")
  pdf(file2write, 
      width = 20,
      height = 16)
  p <- FeaturePlot(object = object2plot, features = markers2plot, 
                   cols = c("grey", "red"), reduction = "umap", label = T)
  print(p)
  dev.off()
  
  # file2write <- paste0(marker_exp_out_path, cell_type_tmp, "_Marker_in_sample_", run_id, "_VlnPlot.pdf")
  # pdf(file2write, 
  #     width = 6,
  #     height = 16)
  # p <- VlnPlot(object2plot, features = markers2plot,
  #              ncol = 1, pt.size = 0.5)
  # print(p)
  # dev.off()
  # 
  # markers2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr == cell_type_tmp]
  # file2write <- paste0(marker_exp_out_path, cell_type_tmp, "_Marker_in_sample_", run_id, "_FeaturePlot_by_sample.pdf")
  # pdf(file2write, 
  #     width = 12,
  #     height = 6)
  # p <- FeaturePlot(object = object2plot, features = markers2plot, 
  #                  cols = c("grey", "red"), reduction = "umap", label = T, split.by = "orig.ident", coord.fixed = T)
  # print(p)
  # dev.off()
}

###########################################
######## Re-Label Cells
###########################################
current.cluster.ids <- c('TAM_C0', "TAM_C2", 'TAM_C6',
                         'CD4Tcell_C8', 'TAM_C9', 'TAM_C10',
                         'Proliferating_TAM_C13',"TAM_C14", "TAM_C15")
new.cluster.ids<- c('MRC1_MacroMono_C0', 
                    "ITGAX_MacroMono_C2", 
                    'MacroMono_C6',
                    'Tcell_C8', 
                    'CDH2_ccRCC_C9', 
                    'PDL1_PDL2_Macrophage_C10',
                    'Proliferating_MacroMono_C13', ## PTPRC/CD45, MRC1/CD206+, CSF1R/CD115+
                    "Macrophage_C14", 
                    "PDL1_PDL2_pDC_C15")
renal.immune.object@meta.data$cell_type <- plyr::mapvalues(as.character(Idents(renal.integrated)), from= current.cluster.ids,to = new.cluster.ids)
Idents(renal.immune.object) <- renal.immune.object@meta.data$cell_type
renal.immune.object@meta.data %>% head()

