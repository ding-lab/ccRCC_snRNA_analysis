# Yige Wu @WashU Feb 2022

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
## input barcode2UMAP
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/recluster/recluster_tumorcells/downsample_fixednumber_and_recluster_tumor_cells_in_selected_samples_katmai/20220222.v1/UMAPData.2000TumorCellReclustered.20220222.v1.tsv")

# plot by each aliquot ----------------------------------------------------
## make different output files
dir_out_png <- paste0(dir_out, "png", "/")
dir.create(dir_out_png)
# dir_out_pdf <- paste0(dir_out, "pdf", "/")
# dir.create(dir_out_pdf)

for (aliquot_tmp in unique(barcode2umap_df$easy_id)) {
  plot_data_df <- barcode2umap_df %>%
    filter(easy_id == aliquot_tmp) %>%
    mutate(Name_TumorCluster = paste0("C", seurat_clusters))
  ## make color for each cluster
  uniq_cluster_colors <- Polychrome::dark.colors(n = length(unique(plot_data_df$Name_TumorCluster)))
  names(uniq_cluster_colors) <- unique(plot_data_df$Name_TumorCluster)
  
  ## make plot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Name_TumorCluster), shape = 16, alpha = 0.8, size = 0.5)
  p <- p + scale_color_manual(values = uniq_cluster_colors, na.translate = T)
  p <- p + theme_bw()
  p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
  p <- p + theme(legend.position = "top")
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank())
  p <- p + theme(axis.line = element_line(colour = "black"),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_text(size = 15))
  p <- p + theme(plot.title = element_text(size = 18))
  p
  ## save plot
  png2write <- paste0(dir_out_png, aliquot_tmp, ".TumorCellOnlyClustering.", ".png")
  png(filename = png2write, width = 900, height = 1000, res = 150)
  print(p)
  dev.off()
  # 
  # file2write <- paste0(dir_out_pdf, aliquot_tmp, ".TumorCellOnlyClustering.", ".pdf")
  # pdf(file2write, width = 5, height = 5.5, useDingbats = F)
  # print(p)
  # dev.off()
}
