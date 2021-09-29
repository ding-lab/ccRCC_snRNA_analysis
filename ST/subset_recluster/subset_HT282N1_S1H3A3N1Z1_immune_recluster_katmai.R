# Yige Wu @WashU Apr 2020

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(ggplot2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input seurat object
srat <- readRDS(file = "../ccRCC_ST/Processed_Data/Seurat/Outputs/TWFU-HT282N1-S1H3A3N1Z1_1Bmn1_1_5.0/TWFU-HT282N1-S1H3A3N1Z1_1Bmn1_1_5.0_processed_multiomic.rds")
print("Finish reading RDS file")
## name of the cell group to recluster
table(srat@meta.data$seurat_clusters)
# 0   1   2   3   4   5   6   7   8   9  10
# 776 771 746 365 318 287 242 238 230 134  79
clusters2process <- c(7)
cat(paste0("Name of the cell clusters to recluster: ", clusters2process, "\n"))
cat("###########################################\n")
resolution_tmp <- 2

# subset ------------------------------------------------------------------
## subset to non-doublet
aliquot_show <- "HT282N1-S1H3A3N1Z1"
## subset to selected clusters
DefaultAssay(srat) <- "RNA"
srat <- subset(srat, idents = clusters2process)
dim(srat)

# before-clustering processing ------------------------------------------------------------
## get variably expressed genes
cat("Start running SCTransform\n")
srat <- SCTransform(srat, vars.to.regress = c("nCount_RNA","percent.mt","S.Score", "G2M.Score"), return.only.var.genes = F)
cat("###########################################\n")
## keep it consistant with individual processing pipeline
cat("Start running RunPCA\n")
srat <- RunPCA(srat, npcs = num_pc, verbose = FALSE)
cat("###########################################\n")
cat("Start running RunUMAP\n")
srat <- RunUMAP(srat,  reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', dims = 1:num_pc)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(srat) <- "ATAC_MACS2"
srat <- RunTFIDF(srat)
srat <- FindTopFeatures(srat, min.cutoff = 'q0')
srat <- RunSVD(srat)
srat <- RunUMAP(srat, reduction = 'lsi', dims = opt$pc_first:opt$pc_num, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


cat("###########################################\n")
cat("Start running FindNeighbors\n")
srat <- FindMultiModalNeighbors(srat,
                                reduction.list = list("pca", "lsi"),
                                dims.list = list(1:opt$pc_num, opt$pc_first:opt$pc_num))
srat <- RunUMAP(srat, nn.name = "weighted.nn",
                reduction.name = "wnn.umap",
                reduction.key = "wnnUMAP_")
srat <- FindClusters(srat, graph.name = "wsnn", algorithm = 3, verbose = T)

## save output
cat("Start saving the reclustered seurat object\n")
file2write <- paste0(dir_out, aliquot_show, ".", "Immune.", "Reclustered.", "Res", resolution_tmp, ".RDS")
saveRDS(object = srat, file = file2write, compress = T)
cat("###########################################\n")

# write dimplot -----------------------------------------------------------
p <- DimPlot(object = srat)
file2write <- paste0(dir_out, "dimplot.png")
png(file2write, width = 800, height = 700, res = 150)
print(p)
dev.off()

# prepare data ------------------------------------------------------------
## input cell type markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/Immune/Cell_state_markers.txt")
DefaultAssay(srat) <- "RNA"
pct_thres <- 15
avgexp_thres <- 0.1
## get the genes within the cell type marker table
genes2plot <-  intersect(gene2celltype_df$Gene, srat@assays$RNA@counts@Dimnames[[1]])
genes2plot <- unique(genes2plot)
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
plotdata_df$gene_cell_type_group <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Gene_set)

p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = features.plot, y = id, color = expvalue_plot, size = pct.exp), shape = 16)
# p <- p +scale_color_gradient2(midpoint=median(plotdata_df$avg.exp, na.rm = T), low="blue", mid="white",
#                               high="red", space ="Lab" )
p <- p + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral")[1:5]), guide = guide_legend(direction = "horizontal"))
p <- p + scale_size_continuous(range = c(0, 8), name="% Expressed", guide = guide_legend(direction = "horizontal"))
# p <- p  + RotatedAxis()
p <- p + facet_grid(.~gene_cell_type_group, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5), 
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 10, angle = 90))
p <- p + theme(legend.position = "bottom")
file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.NotScaled.png")
png(file = file2write, width = 2500, height = 800, res = 150)
print(p)
dev.off()

# plot scaled -------------------------------------------------------------
p <- DotPlot(object = srat, features = genes2plot_filtered, col.min = 0, assay = "RNA")
p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Gene_set)
p <- p + facet_grid(.~gene_cell_type_group, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"),
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 10, angle = 90))
p <- p + theme(legend.position = "bottom")
file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.Scaled.png")
png(file = file2write, width = 3000, height = 1200, res = 150)
print(p)
dev.off()

