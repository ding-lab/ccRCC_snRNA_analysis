# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
## input UMAP info per barcode
barcode2cluster_df <- fread(input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/FindClusters_30_ccRCC_tumorcells_changeresolutions/20220405.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.20220405.v1.tsv", data.table = F)
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# make plot data----------------------------------------------------------
barcode2cluster_df$sample <- mapvalues(x = barcode2cluster_df$orig.ident, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Aliquot.snRNA.WU))
barcode2cluster_df$clusterid_plot <- barcode2cluster_df$integrated_snn_res.1
plot_data_df <- barcode2cluster_df %>%
  group_by(clusterid_plot) %>%
  summarise(number_cells_bycluster = n())
plot_data_df$perc_cells_bycluster = 100*plot_data_df$number_cells_bycluster/sum(plot_data_df$number_cells_bycluster)
plot_data_df$clusterid_plot <- factor(x = plot_data_df$clusterid_plot)
# make colors -------------------------------------------------------------
clusterids <- sort(unique(plot_data_df$clusterid_plot))
colors_bycluster <- Polychrome::palette36.colors(n = (length(clusterids)+1))[-2]
names(colors_bycluster) <- clusterids

# make plots --------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plot_data_df, 
                         mapping = aes(x = clusterid_plot, y = perc_cells_bycluster, fill = clusterid_plot), stat = "identity")
p <- p + scale_fill_manual(values = clusterids)
p <- p + ylab("% in all tumor cells") + xlab("Seurat cluster number")
p <- p + theme_classic()
p <- p + theme(legend.position="none", aspect.ratio=1)
p <- p + theme(axis.title = element_text(size = 15, color = 'black'), axis.text = element_text(size = 12, color = 'black'))
# ## save as pdf
# file2write <- paste0(dir_out, "cellgroup_on_umap.", ".pdf")
# pdf(file = file2write, width = 8, height = 9, useDingbats = F)
# print(p)
# dev.off()
## save as png
file2write <- paste0(dir_out, "perc_cells_bycluster", ".png")
png(filename = file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()

