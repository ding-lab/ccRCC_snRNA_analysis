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
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20201002.v1/31Aliquot.Barcode2CellType.20201002.v1.tsv", data.table = F)
## input srat object
srat <- readRDS(file = "./Resources/snRNA_Processed_Data/Merged_Seurat_Objects/C3N-01200.Tumor_Segments.Merged.20200319.v1.RDS")
## input cell type markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/combine_pt_with_emt_markers/20200911.v1/Kidney_Specific_EMT_Genes.20200911.v1.tsv")
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# specify thresholds ------------------------------------------------------
## filter for genes that are expressed in >25% of one cluster at least
pct_thres <- 20
avgexp_thres <- 1
aliquot_show <- "C3N-01200 Merged "
aliquot2process <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Case == "C3N-01200"]

# map cell type ------------------------------------------------------------------
## subset barcode2celltype
barcode2celltype_filtered_df <- barcode2celltype_df %>%
  filter(orig.ident %in% aliquot2process)
## annotate the easy id
barcode2celltype_filtered_df$id_aliquot_wu <- mapvalues(x = barcode2celltype_filtered_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## identify major cell type in each cluster
count_bycelltype_bycluster_df <- barcode2celltype_filtered_df %>%
  select(id_aliquot_wu, Cell_group7, Cell_group15) %>%
  table() %>%
  data.frame() %>%
  filter(Freq > 0) %>%
  mutate(Id_Cluster_CellType = paste0(id_aliquot_wu, "_", Cell_group15)) %>%
  mutate(Keep = (Freq >= 30 & (Cell_group7 %in% c("Tumor cells", "Transitional cells"))))

barcode2celltype_filtered_df <- barcode2celltype_filtered_df %>%
  mutate(Id_Cluster_CellType = paste0(id_aliquot_wu, "_", Cell_group15)) %>%
  mutate(Id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode))
head(barcode2celltype_filtered_df$Id_aliquot_barcode)
barcode2celltype_filtered_df$Keep <- mapvalues(x = barcode2celltype_filtered_df$Id_Cluster_CellType, from = count_bycelltype_bycluster_df$Id_Cluster_CellType, to = as.vector(count_bycelltype_bycluster_df$Keep))
## subset seurat object
### get the barcodes
metadata_df <- srat@meta.data
metadata_df$integrated_barcode <- rownames(metadata_df)
metadata_df <- metadata_df %>%
  mutate(individual_barcode = str_split_fixed(string = integrated_barcode, pattern = "_", n = 3)[,1]) %>%
  mutate(Id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode))
head(metadata_df$Id_aliquot_barcode)
srat@meta.data$Cell_group <- mapvalues(x = metadata_df$Id_aliquot_barcode, from = barcode2celltype_filtered_df$Id_aliquot_barcode, to = as.vector(barcode2celltype_filtered_df$Id_Cluster_CellType))
unique(srat@meta.data$Cell_group)
Idents(srat) <- "Cell_group"
### subset
srat <- subset(srat, idents = count_bycelltype_bycluster_df$Id_Cluster_CellType[count_bycelltype_bycluster_df$Keep == "TRUE"])
dim(srat)


# prepare data ------------------------------------------------------------
## get the genes within the cell type marker table
genes2plot <-  intersect(gene2celltype_df$Gene, srat@assays$RNA@counts@Dimnames[[1]])
genes2plot <- unique(genes2plot)
## get the pct expressed for each gene in each cluster
p <- DotPlot(object = srat, features = genes2plot, col.min = 0)
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

## add facet
plotdata_df$gene_group2 <- paste0(plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Gene_Group2), 
                                  "\n",
                                  "Markers")
plotdata_df$cell_group <- str_split_fixed(string = plotdata_df$id, pattern = "_", n = 2)[,2]

p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = features.plot, y = id, color = expvalue_plot, size = pct.exp), shape = 16)
p <- p + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral")[1:5]), guide = guide_legend(direction = "horizontal", nrow = 2, byrow = T))
p <- p + scale_size_continuous(range = c(0, 8), name="% Expressed", guide = guide_legend(direction = "horizontal"))
p <- p + facet_grid(cell_group~gene_group2, scales = "free", space = "free", drop = T)
p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10))
p <- p + theme(axis.text.y = element_text( face = "bold", size = 12))
p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5), 
               strip.text.x = element_text(angle = 90, vjust = 0.5, size = 12, face = "bold"), 
               strip.text.y = element_text(angle = 0, vjust = 0.5))
p <- p + theme(axis.title = element_blank())
p <- p + ggtitle(paste0(aliquot_show, " Expression of Cell Type Marker Genes"))
p <- p + labs(colour = "Expression value")
p <- p + theme(legend.position = "bottom")
file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.NotScaled.png")
png(file = file2write, width = 1200, height = 700, res = 150)
print(p)
dev.off()