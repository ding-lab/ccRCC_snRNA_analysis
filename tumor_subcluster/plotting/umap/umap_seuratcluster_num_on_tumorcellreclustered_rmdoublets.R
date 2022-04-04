# Yige Wu @WashU Nov 2020
## for each individual sample tumor cell reclustered, plot UMAP

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
## input barcode2UMAP
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_individual_tumorcellreclustered_on_katmai/20210805.v1/MetaData_TumorCellOnlyReclustered.20210805.v1.tsv")
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)

# plot by each aliquot ----------------------------------------------------
## make different output files
dir_out_png <- paste0(dir_out, "png", "/")
dir.create(dir_out_png)
# dir_out_pdf <- paste0(dir_out, "pdf", "/")
# dir.create(dir_out_pdf)

for (easy_id_tmp in unique(barcode2umap_df$easy_id)) {
  scrublets_df <- barcode2scrublet_df %>%
    filter(Aliquot_WU == easy_id_tmp) %>%
    filter(predicted_doublet)
  if (easy_id_tmp %in% barcode2scrublet_df$Aliquot_WU) {
    plot_data_df <- barcode2umap_df %>%
      filter(easy_id == easy_id_tmp) %>%
      filter(!(barcode_tumorcellreclustered %in% scrublets_df$Barcode)) %>%
      mutate(Name_TumorCluster = paste0("C", (seurat_clusters+1)))
  } else {
    plot_data_df <- barcode2umap_df %>%
      filter(easy_id == easy_id_tmp) %>%
      mutate(Name_TumorCluster = paste0("C", (seurat_clusters+1)))
  }

  ## make color for each cluster
  uniq_cluster_colors <- Polychrome::dark.colors(n = length(unique(plot_data_df$Name_TumorCluster)))
  names(uniq_cluster_colors) <- sort(unique(plot_data_df$Name_TumorCluster))
  
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Name_TumorCluster), shape = 16, alpha = 0.8, size = 1)
  p <- p + scale_color_manual(values = uniq_cluster_colors, na.translate = T)
  p <- p + ggtitle(label = paste0("Tumor-cell subclusters for sample ", easy_id_tmp), subtitle = "Seurat assigned")
  p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = NULL, nrow = 1))
  p <- p + theme_void()
  p <- p + theme(legend.position = "bottom")
  p <- p + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank())
  ## save plot
  png2write <- paste0(dir_out_png, easy_id_tmp,  ".png")
  png(filename = png2write, width = 900, height = 1000, res = 150)
  print(p)
  dev.off()
  # 
  # file2write <- paste0(dir_out_pdf, easy_id_tmp, ".TumorCellOnlyClustering.", ".pdf")
  # pdf(file2write, width = 5, height = 5.5, useDingbats = F)
  # print(p)
  # dev.off()
}
