# Yige Wu @WashU Jan 2020
## Re-cluster each cell type group for the tumor-normal integrated dataset

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

# set case id -------------------------------------------------------------
case_id_tmp <- "C3N-01200"

## set cell type group name to be processed
cell_type_group_tmp <- "Nephron_Epithelium"

# input cell type marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20200131.v1.tsv")

# input reclustered seurat object -----------------------------------------
srat <- readRDS(file = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/recluster_cell_type_groups/20200131.v1/", case_id_tmp, ".Tummor_Normal.Integrated.20200117.v1.", cell_type_group_tmp,".StandardReclustered.20200131.v1.RDS"))

## get cluster marker genes
cluster_markers_tab <- srat@misc$markers

cluster_markers_tab.in.celltype_markers <- cluster_markers_tab %>%
  filter(gene %in% gene2cellType_tab$Gene)

## examine cell type markers by cluster
cluster_tmp <- 3
cluster_markers <-  cluster_markers_tab %>%
  filter(cluster == cluster_tmp) %>%
  filter(gene %in% gene2cellType_tab$Gene)

# dotplot -----------------------------------------------------------------
genes2plot <- gene2cellType_tab$Gene
# genes2plot <- intersect(rownames(seurat_obj@assays$RNA@counts), genes2plot)

p <- DotPlot(object = srat, features = genes2plot, col.min = 0)
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

file2write <- paste0(dir_out,  case_id_tmp, ".Tummor_Normal.Integrated.20200117.v1.", cell_type_group_tmp,".StandardReclustered.", run_id, ".png")
png(file = file2write, width = 4500, height = 1800, res = 150)
print(p)
dev.off()

