# Yige Wu @WashU March 2020
## for each individual sample tumor cell reclustered, plot UMAP

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/aes.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input gene to cell type table
gene2celltype_df <- fread(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/Gene2CellType_Tab.20200220.v1.tsv", data.table = F)

# plot by  ----------------------------------------------------
## input srat object
srat_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
srat <- readRDS(file = srat_path)

## get the genes within the cell type marker table
genes2plot <-  intersect(gene2celltype_df$Gene, srat@assays$RNA@counts@Dimnames[[1]])
genes2plot <- unique(genes2plot)
## get the pct expressed for each gene in each cluster
p <- DotPlot(object = srat, features = genes2plot, col.min = 0)
plot_data <- p$data
## transform the dataframe to matrix to better filter out genes with too low expressin
plot_matrix <- dcast(data = plot_data, formula = features.plot ~ id, value.var = "pct.exp")
plot_matrix %>% head()
## filter for genes that are expressed in >25% of one cluster at least
## replot with the filtered genes plus malignant cell marker genes
malignant_markers <- as.vector(gene2celltype_df$Gene[gene2celltype_df$Cell_Type_Group == "Malignant_Nephron_Epithelium"])
genes2plot_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 25) >= 1, "features.plot"])
genes2plot_filtered <- c(genes2plot_filtered, 
                         as.vector(plot_matrix[(plot_matrix$features.plot %in% malignant_markers) & (rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 5) >= 1), "features.plot"]))
genes2plot_filtered <- unique(genes2plot_filtered)
p <- DotPlot(object = srat, features = genes2plot_filtered, col.min = 0)
p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
p <- p  + RotatedAxis()
p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey50"),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 15, face = "bold"),
               strip.placement = "outside")
## save plot
file2write <- paste0(dir_out, aliquot_tmp, ".TumorCellOnlyClustering.Dotplot.", run_id, ".png")
png(file = file2write, width = 4000, height = 1200, res = 150)
print(p)
dev.off()