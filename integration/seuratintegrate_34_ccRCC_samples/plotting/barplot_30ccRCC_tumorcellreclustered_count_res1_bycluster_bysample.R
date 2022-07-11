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
barcode2cluster_df$clusterid_plot <- paste0("MC",(barcode2cluster_df$integrated_snn_res.1 + 1))

plot_data_df <- barcode2cluster_df %>%
  group_by(sample, clusterid_plot) %>%
  summarise(number_cells_bycluster_bysample = n())
plot_data_df2 <- barcode2cluster_df %>%
  group_by(sample) %>%
  summarise(number_cells_bysample = n())
plot_data_df3 <- barcode2cluster_df %>%
  group_by(clusterid_plot) %>%
  summarise(number_cells_bycluster = n())
plot_data_df$number_cells_bysample <- mapvalues(x = plot_data_df$sample, from = plot_data_df2$sample, to = as.vector(plot_data_df2$number_cells_bysample))
plot_data_df$number_cells_bysample <- as.numeric(plot_data_df$number_cells_bysample)
plot_data_df$number_cells_bycluster <- mapvalues(x = plot_data_df$clusterid_plot, from = plot_data_df3$clusterid_plot, to = as.vector(plot_data_df3$number_cells_bycluster))
plot_data_df$number_cells_bycluster <- as.numeric(plot_data_df$number_cells_bycluster)
plot_data_df <- plot_data_df %>%
  mutate(perc_cells_bysample_eachcluster = (number_cells_bycluster_bysample/number_cells_bycluster)*100)
plot_data_df$clusterid_plot <- factor(x = plot_data_df$clusterid_plot, levels = paste0("MC", 1:18))

# make colors -------------------------------------------------------------
sampleids_ordered <- sort(unique(plot_data_df$sample))
colors_cellgroup <- Polychrome::palette36.colors(n = length(sampleids_ordered))
names(colors_cellgroup) <- sampleids_ordered

# make plots --------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plot_data_df, 
                         mapping = aes(x = clusterid_plot, y = perc_cells_bysample_eachcluster, fill = sample), stat = "identity")
p <- p + scale_fill_manual(values = colors_cellgroup)
p <- p + guides(fill = guide_legend(override.aes = list(size=4), title = NULL))
p <- p + ylab("% cells by sample")
p <- p + theme_classic()
p <- p + theme(#axis.ticks.x=element_blank(), axis.line = element_blank(), 
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15, color = "black"),
               axis.text.y = element_text(size = 15, color = "black"), axis.title.x = element_blank(), axis.title.y = element_text(size = 15))
p <- p + theme(legend.position="right", aspect.ratio=1, legend.text = element_text(size = 14))
## save as pdf
file2write <- paste0(dir_out, "perc_cells_bysample_eachcluster", ".pdf")
pdf(file = file2write, width = 8, height = 5, useDingbats = F)
print(p)
dev.off()
## save as png
file2write <- paste0(dir_out, "perc_cells_bysample_eachcluster", ".png")
png(filename = file2write, width = 1200, height = 800, res = 150)
print(p)
dev.off()

#