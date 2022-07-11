# Yige Wu @WashU Dec 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
library(ggrastr)
library(ggplot2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode with UMAP coordinates and manual group id
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id_PC50/20220307.v1/Barcode2TumorSubclusterId.20220307.v1.tsv")

# plot by each aliquot ----------------------------------------------------
## make different output files
dir_out_png <- paste0(dir_out, "png", "/")
dir.create(dir_out_png)
# dir_out_pdf <- paste0(dir_out, "pdf", "/")
# dir.create(dir_out_pdf)
colors_all <- Polychrome::dark.colors(n = 24)

# for (easy_id_tmp in "C3L-00079-T1") {
for (easy_id_tmp in unique(barcode2umap_df$easy_id)) {
  
  plot_data_df <- barcode2umap_df %>%
    filter(easy_id == easy_id_tmp) %>%
    mutate(Name_TumorCluster = str_split_fixed(string = Cluster_Name.cutoff50cells, pattern = "_", n = 2)[,2]) %>%
    mutate(Name_TumorCluster = ifelse(Name_TumorCluster == "CNA", "Minor cluster (<50 cells)", Name_TumorCluster))

  ## make color for each cluster
  names_cluster_tmp <- unique(plot_data_df$Name_TumorCluster)
  length_clusters <- length(names_cluster_tmp)
  if("Minor cluster (<50 cells)" %in% names_cluster_tmp) {
    uniq_cluster_colors <- c(colors_all[1:(length_clusters-1)], "grey40")
  } else {
    uniq_cluster_colors <-colors_all[1:length_clusters]
  }
  names(uniq_cluster_colors) <- names_cluster_tmp
  
  ## make plot
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Name_TumorCluster), shape = 16, alpha = 0.8, size = 1)
  p <- p + scale_color_manual(values = uniq_cluster_colors, na.translate = T)
  # p <- p + ggtitle(label = easy_id_tmp, subtitle = "Manually grouped tumor clusters")
  p <- p + ggtitle(label = easy_id_tmp, subtitle = "50 PCs")
  p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = NULL, nrow = 2, label.theme = element_text(size = 20)))
  p <- p + theme_void()
  p <- p + theme(legend.position = "bottom")
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), title = element_text(size = 25))
  # save plot
  png2write <- paste0(dir_out_png, easy_id_tmp, ".png")
  png(filename = png2write, width = 900, height = 1100, res = 150)
  print(p)
  dev.off()
  
  # file2write <- paste0(dir_out_pdf, easy_id_tmp,".pdf")
  # pdf(file2write, width = 5, height = 5.5, useDingbats = F)
  # print(p)
  # dev.off()
}
