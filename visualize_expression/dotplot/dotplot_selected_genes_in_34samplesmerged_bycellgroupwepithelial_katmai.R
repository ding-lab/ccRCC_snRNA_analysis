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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/merging/merge_34_ccRCC_samples/20211005.v1//ccRCC.34samples.Merged.20211005.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")
## spcify assay
assay_process <- "SCT"
slot_process <- "data"
cat(paste0("Assay: ", assay_process, "\n"))
pct_thres <- 15
avgexp_thres <- 0.1
## specify genes to plot
genes_plot <- c("CA9", "PAX8", "PAX2", "CD24", "LRP2")
# set ident ---------------------------------------------------------------
## make unique id for each barcode in the cell type table
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(barcode = paste0(orig.ident, "_", individual_barcode)) %>%
  mutate(cell_group = Cell_group_w_epithelialcelltypes)
## 
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
## remove unknown cell types
cellgroups_keep <- unique(barcode2celltype_df$cell_group); cellgroups_keep <- cellgroups_keep[!(cellgroups_keep %in% c("Unknown"))]
srat <- subset(x = srat, idents = cellgroups_keep)
dim(srat)

# plot not scaled -------------------------------------------------------------
## get the pct expressed for each gene in each cluster
# p <- DotPlot(object = srat, features = genes_plot, col.min = 0, assay = "RNA")
# expdata_df <- p$data
# plotdata_df <- expdata_df %>%
#   filter(features.plot %in% genes_plot)
# expvalue_top <- quantile(x = plotdata_df$avg.exp, probs = 0.95)
# plotdata_df <- plotdata_df %>%
#   mutate(expvalue_plot = ifelse(avg.exp >= expvalue_top, expvalue_top, avg.exp))
# p <- ggplot()
# p <- p + geom_point(data = plotdata_df, mapping = aes(x = features.plot, y = id, color = expvalue_plot, size = pct.exp), shape = 16)
# # p <- p +scale_color_gradient2(midpoint=median(plotdata_df$avg.exp, na.rm = T), low="blue", mid="white",
# #                               high="red", space ="Lab" )
# p <- p + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral")[1:5]), guide = guide_legend(direction = "horizontal"))
# p <- p + scale_size_continuous(range = c(0, 8), name="% Expressed", guide = guide_legend(direction = "horizontal"))
# # p <- p  + RotatedAxis()
# p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
#                panel.border = element_rect(color = "black", fill = NA, size = 0.5),
#                panel.background = element_blank())
# p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5), 
#                strip.text.x = element_text(angle = 0, vjust = 0.5),
#                strip.text.y = element_text(angle = 0, vjust = 0.5),
#                axis.text.x = element_text(size = 10, angle=90,hjust=0.95,vjust=0.2))
# p <- p + theme(legend.position = "bottom")
# file2write <- paste0(dir_out, "CellTypeMarkerExp.NotScaled.png")
# png(file = file2write, width = 1000, height = 1000, res = 150)
# print(p)
# dev.off()

# plot scaled -------------------------------------------------------------
p <- DotPlot(object = srat, features = genes_plot, col.min = 0, assay = "RNA")
p <- p + theme(axis.text.x = element_text(size = 10, angle=90,hjust=0.95,vjust=0.2), axis.title = element_blank())
# p <- p + theme(legend.position = "bottom")
p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
               panel.border = element_rect(color = "black", fill = NA, size = 0.5),
               panel.background = element_blank())
p <- p + theme(strip.background = element_rect(color = NA, fill = NA, size = 0.5), 
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 10, angle=90,hjust=0.95,vjust=0.2))
file2write <- paste0(dir_out, "CellTypeMarkerExp.Scaled.png")
png(file = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()


