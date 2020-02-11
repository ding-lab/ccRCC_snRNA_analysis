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
case_id_tmp <- "C3N-01200"
  
# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/scripts/src/Marker_Gene_Tables/Gene2CellType_Tab.20200203.v1.tsv", data.table = F)

# plot dotplot by sample --------------------------------------------------
seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_tumor_normal/20200117.v1/", case_id_tmp, ".Tummor_Normal.Integrated.20200117.v1.RDS")
deg_tab_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/run_deg_same_patient_tumor_normal/20200117.v1/", case_id_tmp, ".Tumor_Normal.DEGs.Pos.txt")

## input seurat object
seurat_obj <- readRDS(file = seurat_obj_path)
DefaultAssay(seurat_obj) <- "RNA"

## plot all marker gene expression
genes2plot <- gene2cellType_tab$Gene
## get the pct expressed for each gene in each cluster
p <- DotPlot(object = seurat_obj, features = genes2plot, col.min = 0)
plot_data <- p$data
## transform the dataframe to matrix to better filter out genes with too low expressin
library(data.table)
plot_matrix <- dcast(data = plot_data, formula = features.plot ~ id, value.var = "pct.exp")
plot_matrix %>% head()
## filter for genes that are expressed in >25% of one cluster at least
plot_matrix_filtered <- plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 25) >= 1,]
## replot with the filtered genes
p <- DotPlot(object = seurat_obj, features = plot_matrix_filtered$features.plot, col.min = 0)
p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Group)
p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type1)
p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type2)
p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type3)
p <- p  + RotatedAxis()
p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey50"),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 15, face = "bold"),
               strip.placement = "outside")

file2write <- paste0(dir_out,  case_id_tmp,".Tumor_Normal.Dotplot.CellTypeMarkers.", run_id, ".png")
png(file = file2write, width = 5500, height = 1800, res = 150)
print(p)
dev.off()



