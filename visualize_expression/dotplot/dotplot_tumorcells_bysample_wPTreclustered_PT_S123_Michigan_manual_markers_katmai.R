# Yige Wu @WashU Aug 2020
## finding differentially expressed gene for each cell type using integrared object

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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/merging/merge_35_samples/20210802.v1/RCC.35samples.Merged.20210802.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_byepithelialcelltypes_bysample_wPTreclustered/20210809.v1/Barcode_byEpithelialCelltypes_BySample.20210809.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")
## input cell type markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Human.Gene2CellType.20210819.tsv")

## spcify assay
assay_process <- "SCT"
slot_process <- "data"
cat(paste0("Assay: ", assay_process, "\n"))
pct_thres <- 0
avgexp_thres <- 0

# set ident ---------------------------------------------------------------
srat@meta.data$individual_barcode <- str_split_fixed(string = rownames(srat@meta.data), pattern = "_", n = 2)[,1]
## check if the individual_barcode is mapped right
srat@meta.data %>% head()
## add cell id to the seurat meta data
srat@meta.data$id_cell <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$individual_barcode)
srat@meta.data$cell_group_process <- mapvalues(x = srat@meta.data$id_cell, from = barcode2celltype_df$barcode, to = as.vector(barcode2celltype_df$cell_group))
unique(srat@meta.data$cell_group_process)
Idents(srat) <- "cell_group_process" 
dim(srat)

# subset object -----------------------------------------------------------
count_bycellgroup_df <- barcode2celltype_df %>%
  select(cell_group) %>%
  table() %>%
  as.data.frame() %>%
  rename(cell_group = ".") %>%
  mutate(cell_type = str_split_fixed(string = cell_group, pattern = "_", n = 2)[,1])
count_bycellgroup_keep_df <- count_bycellgroup_df %>%
  filter(Freq >= 50) %>%
  filter(cell_type %in% c("EMT tumor cells", "PT", "Tumor cells")) %>%
  filter(!grepl(x = cell_group, pattern = "Unknown"))
cellgroups_keep <- count_bycellgroup_keep_df$cell_group; cellgroups_keep <- as.vector(cellgroups_keep)
srat <- subset(x = srat, idents = cellgroups_keep)
dim(srat)

# prepare data ------------------------------------------------------------
gene2celltype_df <- gene2celltype_df %>%
  filter(Gene %in% c("LRP2", "SLC5A12", "SLC5A2", "ACSM3", "CA4", "GATM", "SLC3A1", "GLYAT", "ITGB8", "ALPK2", "CFH", "KLK6", "KRT19", "ALDOB", "PDZK1IP1", "SLC7A13")) %>%
  select(Gene, Cell_Type2) %>%
  mutate(Cell_Type2 = gsub(pattern = "PT S", replacement = "S", x = Cell_Type2))
gene2celltype_df$Cell_Type2[gene2celltype_df$Gene %in% c("KRT19", "LRP2", "ALDOB", "PDZK1IP1")] <- "PT"
gene2celltype_df$Cell_Type2[gene2celltype_df$Gene %in% c("GLYAT")] <- "PT-A"
genes2celltype_add_df <- data.frame(Gene = c("CYP24A1", "HKDC1", "SLC26A3", "NDC80"))
genes2celltype_add_df$Cell_Type2 <- "Other"
gene2celltype_df <- rbind(gene2celltype_df, genes2celltype_add_df)

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
plotdata_df$gene_cell_type2 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
plotdata_df$gene_cell_type2 <- factor(x = plotdata_df$gene_cell_type2, levels = c("PT", "S1", "S1/S2", "S2", "S3", "PT-A", "PT-B", "PT-C", "Other"))

plotdata_df$cell_type <- plyr::mapvalues(x = plotdata_df$id, from = count_bycellgroup_keep_df$cell_group, to = as.vector(count_bycellgroup_keep_df$cell_type))
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = features.plot, y = id, color = expvalue_plot, size = pct.exp), shape = 16)
# p <- p +scale_color_gradient2(midpoint=median(plotdata_df$avg.exp, na.rm = T), low="blue", mid="white",
#                               high="red", space ="Lab" )
p <- p + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral")[1:5]), guide = guide_legend(direction = "horizontal"))
p <- p + scale_size_continuous(range = c(0, 8), name="% Expressed", guide = guide_legend(direction = "horizontal"))
# p <- p  + RotatedAxis()
# p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + facet_grid(cell_type~gene_cell_type2, scales = "free", space = "free", drop = T)

p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5), 
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 10, angle=90,hjust=0.95,vjust=0.2))
p <- p + theme(legend.position = "bottom")
file2write <- paste0(dir_out, "CellTypeMarkerExp.NotScaled.png")
png(file = file2write, width = 1300, height = 2000, res = 150)
# png(file = file2write, width = 1200, height = 2000, res = 150)
print(p)
dev.off()

# plot scaled -------------------------------------------------------------
p <- DotPlot(object = srat, features = genes2plot_filtered, col.min = 0, assay = "RNA")
p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
p$data$cell_type <- plyr::mapvalues(x = p$data$id, from = count_bycellgroup_keep_df$cell_group, to = as.vector(count_bycellgroup_keep_df$cell_type))
p$data$gene_cell_type2 <- factor(x = p$data$gene_cell_type2, levels = c("PT", "S1", "S1/S2", "S2", "S3", "PT-A", "PT-B", "PT-C", "Other"))
# p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + facet_grid(cell_type~gene_cell_type2, scales = "free", space = "free", drop = T)

p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"),
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 10, angle=90,hjust=0.95,vjust=0.2))
p <- p + theme(legend.position = "bottom")
file2write <- paste0(dir_out, "CellTypeMarkerExp.Scaled.png")
png(file = file2write, width = 1500, height = 2000, res = 150)
print(p)
dev.off()


