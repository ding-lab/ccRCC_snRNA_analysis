# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
integrated_umap_df <- fread(input = "./Resources/Analysis_Results/fetch_data/fetchdata_34_ccRCC_samples_merged_katmai/20211005.v1/ccRCC.34Sample.Merged.Metadata.20211005.v1.tsv", data.table = F)
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# make plot data----------------------------------------------------------
plot_data_df <- integrated_umap_df
plot_data_df$sample <- mapvalues(x = plot_data_df$orig.ident, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Aliquot.snRNA.WU))

# make colors -------------------------------------------------------------
sampleids_ordered <- unique(plot_data_df$sample)
colors_cellgroup <- Polychrome::palette36.colors(n = length(sampleids_ordered))
# swatch(colors_cellgroup)
# swatch(Polychrome::palette36.colors(n = 36))
names(colors_cellgroup) <- sampleids_ordered

# make plots --------------------------------------------------------------
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = sample),
                         alpha = 1, size = 0.1, shape = 16)
p <- p + scale_color_manual(values = colors_cellgroup)
p <- p + guides(colour = guide_legend(override.aes = list(size=4), title = NULL))
p <- p + theme_void()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_blank())
               # axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position="right", aspect.ratio=1)
p
# ## save as pdf
# file2write <- paste0(dir_out, "cellgroup_on_umap.", ".pdf")
# pdf(file = file2write, width = 8, height = 9, useDingbats = F)
# print(p)
# dev.off()
## save as png
file2write <- paste0(dir_out, "cellgroup_on_umap.", ".png")
png(filename = file2write, width = 1500, height = 1000, res = 150)
print(p)
dev.off()

#