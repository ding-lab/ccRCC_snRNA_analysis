# Yige Wu @WashU Apr 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input seurat object
srat <- readRDS(file = "./Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0012280004/pf3000_fmin200_fmax10000_cmin3000_cmax10000_mito_max0.1/CPT0012280004_processed.rds")
## input cell type markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Human.Gene2CellType.20210223.tsv")
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/individual_sample/map_celltype_for_C3N-00317-T1/20210305.v1/C3N-00317-T1.Barcode2CellType.20210305.v1.tsv", data.table = F)

# specify thresholds ------------------------------------------------------
## filter for genes that are expressed in >25% of one cluster at least
aliquot_show <- "C3N-00317-T1"
pct_thres <- 20
avgexp_thres <- 0.1

# prepare data ------------------------------------------------------------
## get the genes within the cell type marker table
genes2plot <-  intersect(gene2celltype_df$Gene, srat@assays$RNA@counts@Dimnames[[1]])
genes2plot <- unique(genes2plot)
## count cell group cell number
BC <- srat@meta.data %>% rownames
srat@meta.data$Cell_type.detailed <- mapvalues(x = BC, from = barcode2celltype_df$individual_barcode, to = as.vector(barcode2celltype_df$Cell_type.detailed))
srat@meta.data$Cell_group <- paste0(srat@meta.data$seurat_clusters, "_", srat@meta.data$Cell_type.detailed)
count_bycelltype_bycluster_df <- srat@meta.data %>%
  as.data.frame() %>%
  select(Cell_group) %>%
  table() %>%
  data.frame() %>%
  filter(Freq > 0) %>%
  mutate(Keep = (Freq >= 30)) %>%
  rename(Cell_group = '.')
## change cell idents
Idents(srat) <- "Cell_group"
## subset seurat object
srat <- subset(srat, idents = count_bycelltype_bycluster_df$Cell_group[count_bycelltype_bycluster_df$Keep == "TRUE"])

## get the pct expressed for each gene in each cluster
p <- DotPlot(object = srat, features = genes2plot, col.min = 0, assay = "RNA")
expdata_df <- p$data
## filter genes based on the percentage expressed
pct_matrix <- dcast(data = expdata_df, formula = features.plot ~ id, value.var = "pct.exp")
genes_pct_filtered <- as.vector(pct_matrix[rowSums(pct_matrix[,unique(as.vector(expdata_df$id))] > pct_thres) >= 1, "features.plot"])
## filter genes based on the average expression
avgexp_matrix <- dcast(data = expdata_df, formula = features.plot ~ id, value.var = "avg.exp")
genes_exp_filtered <- as.vector(avgexp_matrix[rowSums(avgexp_matrix[,unique(as.vector(expdata_df$id))] > avgexp_thres) >= 1, "features.plot"])
## intersect
genes2plot_filtered <- intersect(unique(genes_exp_filtered), unique(genes_pct_filtered))

# plot not scaled -------------------------------------------------------------
plotdata_df <- expdata_df %>%
  filter(features.plot %in% genes2plot_filtered)
expvalue_top <- quantile(x = plotdata_df$avg.exp, probs = 0.95)
plotdata_df <- plotdata_df %>%
  mutate(expvalue_plot = ifelse(avg.exp >= expvalue_top, expvalue_top, avg.exp))
plotdata_df$gene_cell_type_group <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
plotdata_df$gene_cell_type1 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
plotdata_df$gene_cell_type2 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
plotdata_df$gene_cell_type3 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
plotdata_df$gene_cell_type4 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)

p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = features.plot, y = id, color = expvalue_plot, size = pct.exp), shape = 16)
# p <- p +scale_color_gradient2(midpoint=median(plotdata_df$avg.exp, na.rm = T), low="blue", mid="white",
#                               high="red", space ="Lab" )
p <- p + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral")[1:5]), guide = guide_legend(direction = "horizontal"))
p <- p + scale_size_continuous(range = c(0, 8), name="% Expressed", guide = guide_legend(direction = "horizontal"))
# p <- p  + RotatedAxis()
p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5), 
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 10, angle = 90))
p <- p + theme(legend.position = "bottom")
file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.NotScaled.png")
png(file = file2write, width = 3000, height = 1500, res = 150)
print(p)
dev.off()

# plot scaled -------------------------------------------------------------
p <- DotPlot(object = srat, features = genes2plot_filtered, col.min = 0, assay = "RNA")
p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"),
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 10, angle = 90))
p <- p + theme(legend.position = "bottom")
file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.Scaled.png")
png(file = file2write, width = 2500, height = 800, res = 150)
print(p)
dev.off()
