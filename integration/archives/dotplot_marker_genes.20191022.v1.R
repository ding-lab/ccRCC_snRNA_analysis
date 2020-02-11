# Yige Wu @WashU Sep 2019
## for integrating two snRNA datasets for sample CPT0086820004 and CPT0075130004 (from cellranger output with premrna reference)

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

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191021.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191022.v4.tsv")

# input DEG ---------------------------------------------------------------
renal.markers <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/integration/integrate_seurat_objects/20191015.v1/Renal.DEGs.Pos.txt", data.table = F)

# plot DEG in dot plot before cell type assignment-----------------------------------------------------------
object2plot <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/integration/integrate_seurat_objects/20191021.v1/Renal_Integrated.20191021.v1.RDS")
DefaultAssay(object2plot) <- "RNA"

genes2plot <- renal.markers$gene %>% unique()
genes2plot <- genes2plot[genes2plot %in% gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]]
genes2plot <- c("CA9", "HIF1A", "EPAS1", "VHL", "VIM1", "MKI67", "CASP3", genes2plot)
genes2plot

p <- DotPlot(object = object2plot, features = genes2plot, col.min = 0, dot.min = 0.00000001)
p$data$cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
p$data$Is_Immune_Cell_Type_Marker <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Is_Immune_Cell_Type_Marker)
p <- p  + RotatedAxis()
p <- p + facet_grid(~Is_Immune_Cell_Type_Marker + cell_type, scales = "free_x", space = "free", shrink = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.placement = "outside")
file2write <- paste0(dir_out, "Dotplot_AllClusters_DEG_Marker_Exp_without_CellType.", run_id, ".pdf")
pdf(file2write, width = 20, height = 9)
print(p)
dev.off()

file2write <- paste0(dir_out, "Dotplot_AllClusters_DEG_Marker_Exp_without_CellType.", run_id, ".png")
png(file = file2write, width = 3000, height = 1500, res = 150)
print(p)
dev.off()

# plot all markers in dot plot before cell type assignment-----------------------------------------------------------
DefaultAssay(object2plot) <- "RNA"

genes2plot <- gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]
genes2plot <- c("CA9", "HIF1A", "EPAS1", "VHL", "VIM1", "MKI67", "CASP3", genes2plot)
genes2plot <- intersect(genes2plot, rownames(object2plot@assays$RNA@counts))
genes2plot

p <- DotPlot(object = object2plot, features = genes2plot, col.min = 0, dot.min = 0.00000001)
p$data$cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
p <- p  + RotatedAxis()
p <- p + facet_grid(.~cell_type, scales = "free_x", space = "free", shrink = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.placement = "outside")
file2write <- paste0(dir_out, "Dotplot_AllClusters_Marker_Exp_without_CellType.", run_id, ".pdf")
pdf(file2write, width = 35, height = 7)
print(p)
dev.off()



# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191022.v1.tsv", data.table = F)

# plot dot plot after cell type assignment-----------------------------------------------------------
nonimmune_clusters <- cluster2celltype_tab %>%
  filter(Is_Immune == "No") %>%
  select(Cluster)
nonimmune_clusters <- nonimmune_clusters$Cluster

genes2plot <- renal.markers$gene %>% unique()
genes2plot <- genes2plot[genes2plot %in% gene2cellType_tab$Gene[gene2cellType_tab$Cell_Type_Abbr != "Other"]]
genes2plot <- c("CA9", "MKI67", "VIM", "CASP3", genes2plot)
genes2plot

p <- DotPlot(object = object2plot, features = genes2plot, col.min = 0)
p$data$gene_cell_type <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Abbr)
p$data$cluster_cell_type <- plyr::mapvalues(p$data$id, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
p$data$gene_cell_type[p$data$features.plot == "CASP3"] <- "Other"
p$data$is_immune <- ifelse(p$data$id %in% nonimmune_clusters, "Non_Immune", "Immune")
p$data$is_immune <- factor(p$data$is_immune, levels = c("Non_Immune", "Immune"))
p$data$Is_Immune_Cell_Type_Marker <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Is_Immune_Cell_Type_Marker)

p <- p  + RotatedAxis()
p <- p + facet_grid(is_immune + cluster_cell_type~Is_Immune_Cell_Type_Marker+gene_cell_type, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 9),
               strip.placement = "outside")

file2write <- paste0(dir_out, "Dotplot_AllClusters_Marker_Exp.", run_id, ".png")
png(file = file2write, width = 3000, height = 1500, res = 150)
print(p)
dev.off()
