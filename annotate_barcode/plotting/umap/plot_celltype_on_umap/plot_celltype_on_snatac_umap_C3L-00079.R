# Yige Wu @WashU Apr 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggrastr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input UMAP info per barcode
umap_df <- fread(input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Barcode_Annotation/UMAP/UMAP_data_C3L_00079_T1.20201113.tsv", data.table = F)
aliquot_show <- "C3L-00079"

# plot for cell group----------------------------------------------------------
plotdata_df <- umap_df
plotdata_df <- plotdata_df %>%
  dplyr::rename(Cell_group = cell_group_manual)

p <- ggplot()
p <- p + geom_point_rast(data = plotdata_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                         alpha = 1, size = 0.5)
p <- p + scale_color_manual(values = cellgroup_colors[unique(plotdata_df$Cell_group)])
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(axis.line=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank())
p <- p + theme(axis.title = element_text(size = 20))
p <- p + theme(legend.position = "none")
p
file2write <- paste0(dir_out, aliquot_show, ".nolegend.pdf")
pdf(file2write, width = 8, height = 8, useDingbats = F)
print(p)
dev.off()

p <- ggplot()
p <- p + geom_point_rast(data = plotdata_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                         alpha = 1, size = 0.5)
p <- p + scale_color_manual(values = cellgroup_colors[unique(plotdata_df$Cell_group)])
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(axis.line=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank())
p <- p + theme(axis.title = element_text(size = 20))
p <- p + guides(colour = guide_legend(override.aes = list(size=8), nrow = 1, byrow = T))
p <- p + theme(legend.position = "bottom", legend.text = element_text(size = 20), legend.title = element_text(size = 25))
p
file2write <- paste0(dir_out, aliquot_show, ".withlegend.pdf")
pdf(file2write, width = 10, height = 9, useDingbats = F)
print(p)
dev.off()
