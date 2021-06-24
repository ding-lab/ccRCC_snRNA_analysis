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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type per barcode table
# barcode2celltype_df <- fread(input = "./Data_Freezes/V1/snRNA/Cell_Type_Assignment/31Aliquot.Barcode2CellType.20201027.v1.tsv", data.table = F)
barcode2celltype_df <- fread(input = "./Data_Freezes/V2/snRNA/Cell_Type_Assignment/33Aliquot.Barcode2CellType.20210423.v1.tsv", data.table = F)
## input srat object
srat <- readRDS(file = "./Data_Freezes/V1/snRNA/Merged_Seurat_Objects/C3N-01200.Tumor_Segments.Merged.20200319.v1.RDS")
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input cell type markers
gene2celltype_df <- data.frame(Gene = c("MRC1", "CD163", "MSR1",
                                        "CD86", "FCGR3A", "TLR2", "SRSF10", "SIGLEC1"),
                               Gene_Group = c(rep(x = "M2", 3), rep(x = "M1", 5)))

# specify thresholds ------------------------------------------------------
## filter for genes that are expressed in >25% of one cluster at least
pct_thres <- 20
avgexp_thres <- 1
aliquot_show <- "C3N-01200 Merged"
aliquot2process <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Case == "C3N-01200"]

# map cell type ------------------------------------------------------------------
## subset barcode2celltype
barcode2celltype_filtered_df <- barcode2celltype_df %>%
  filter(orig.ident %in% aliquot2process) %>%
  mutate(Id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode))
## annotate the easy id
barcode2celltype_filtered_df$id_aliquot_wu <- mapvalues(x = barcode2celltype_filtered_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## subset seurat object
### get the barcodes
metadata_df <- srat@meta.data
metadata_df$integrated_barcode <- rownames(metadata_df)
metadata_df <- metadata_df %>%
  mutate(individual_barcode = str_split_fixed(string = integrated_barcode, pattern = "_", n = 3)[,1]) %>%
  mutate(Id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode))
head(metadata_df$Id_aliquot_barcode)
srat@meta.data$Cell_group <- mapvalues(x = metadata_df$Id_aliquot_barcode, from = barcode2celltype_filtered_df$Id_aliquot_barcode, to = as.vector(barcode2celltype_filtered_df$Cell_group13))
srat@meta.data$id_aliquot_wu <- mapvalues(x = srat@meta.data$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
unique(srat@meta.data$Cell_group)
Idents(srat) <- "Cell_group"
### subset
srat <- subset(srat, idents = "Macrophages")
dim(srat)
Idents(srat) <- "id_aliquot_wu"
srat <- subset(srat, idents = c("C3N-01200-T1","C3N-01200-T2","C3N-01200-T3"))

# prepare data ------------------------------------------------------------
## get the genes within the cell type marker table
genes2plot <-  intersect(gene2celltype_df$Gene, srat@assays$RNA@counts@Dimnames[[1]])
genes2plot <- unique(genes2plot); genes2plot
## get the pct expressed for each gene in each cluster
p <- DotPlot(object = srat, features = genes2plot, col.min = 0)
expdata_df <- p$data

# plot scaled -------------------------------------------------------------
## add facet
p$data$gene_group <- paste0(plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = as.vector(gene2celltype_df$Gene_Group)), 
                                 " Markers")
p$data$id <- factor(x = as.vector(p$data$id), levels =  rev(c("C3N-01200-T1","C3N-01200-T2","C3N-01200-T3")))
p <- p + facet_grid(.~gene_group, scales = "free", space = "free", drop = T)
p <- p + theme(axis.text.x = element_text(angle = 90, face = "italic", size = 12))
p <- p + theme(axis.text.y = element_text(size = 15))
p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5), 
               strip.text.x = element_text(angle = 0, vjust = 0.5, size = 12), 
               strip.text.y = element_text(angle = 0, vjust = 0.5))
p <- p + theme(axis.title = element_blank())
p <- p + guides(size = guide_legend(title.position = "top", 
                                    nrow = 1, label.theme = element_text(size = 14)))
# p <- p + ggtitle(paste0(aliquot_show, " Macrophages"))
p <- p + labs(colour = "Expression value")
p <- p + theme(legend.position = "bottom")
# file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.Scaled.png")
# png(file = file2write, width = 800, height = 500, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.Scaled.pdf")
# pdf(file2write, width = 4, height = 2.5, useDingbats = F)
pdf(file2write, width = 4, height = 2.5, useDingbats = F)
print(p)
dev.off()
# file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.Scaled.legend.pdf")
# pdf(file2write, width = 8, height = 4, useDingbats = F)
# print(p)
# dev.off()

# # plot not scaled -------------------------------------------------------------
# plotdata_df <- expdata_df
# expvalue_top <- quantile(x = plotdata_df$avg.exp, probs = 0.95)
# plotdata_df <- plotdata_df %>%
#   mutate(expvalue_plot = ifelse(avg.exp >= expvalue_top, expvalue_top, avg.exp))
# 
# ## add facet
# plotdata_df$gene_group <- paste0(plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = as.vector(gene2celltype_df$Gene_Group)), 
#                                   " Markers")
# plotdata_df$id <- factor(x = as.vector(plotdata_df$id), levels =rev(c("C3N-01200-T1","C3N-01200-T2","C3N-01200-T3")))
# 
# p <- ggplot()
# p <- p + geom_point(data = plotdata_df, mapping = aes(x = features.plot, y = id, color = expvalue_plot, size = pct.exp), shape = 16)
# p <- p + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral")[1:5]), guide = guide_legend(direction = "horizontal", nrow = 2, byrow = T))
# p <- p + scale_size_continuous(range = c(0, 8), name="% Expressed", guide = guide_legend(direction = "horizontal"))
# p <- p + facet_grid(.~gene_group, scales = "free", space = "free", drop = T)
# p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10))
# p <- p + theme(axis.text.y = element_text( face = "bold", size = 12))
# p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
#                panel.border = element_rect(color = "black", fill = NA, size = 0.5),
#                panel.background = element_blank())
# p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5), 
#                strip.text.x = element_text(angle = 0, vjust = 0.5, size = 12, face = "bold"), 
#                strip.text.y = element_text(angle = 0, vjust = 0.5))
# p <- p + theme(axis.title = element_blank())
# # p <- p + ggtitle(paste0(aliquot_show, " Macrophages"))
# p <- p + labs(colour = "Expression value")
# p <- p + theme(legend.position = "bottom")
# file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.NotScaled.png")
# png(file = file2write, width = 800, height = 500, res = 150)
# print(p)
# dev.off()