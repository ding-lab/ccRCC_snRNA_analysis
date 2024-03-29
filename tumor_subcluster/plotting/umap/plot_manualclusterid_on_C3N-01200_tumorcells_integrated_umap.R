# Yige Wu @WashU Mar 2022

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
  "Seurat",
  "ggrastr",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode2UMAP
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/integrate_C3N-01200_tumorcells/integrate_C3N-01200_tumorcells/20220608.v1/C3N-01200.Tumorcells.Integrated.UMAP_data.20220608.v1.tsv")
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)

# pre-process --------------------------------------------------------------
barcode2umap_df <- merge(x = barcode2umap_df %>%
                           mutate(barcode = str_split_fixed(string = barcode_merged, pattern = "_", n = 2)[,1]), 
                         y = barcode2tumorsubcluster_df, 
                         by.x = c("orig.ident", "barcode"),
                         by.y = c("orig.ident", "barcode"),
                         all.x = T)

# plot----------------------------------------------------
plot_data_df <- barcode2umap_df %>%
  mutate(Text_Cluster = gsub(x = Cluster_Name, pattern = "C3N\\-01200\\-", replacement = ""))
## make colors
texts_cluster_uniq <- unique(plot_data_df$Text_Cluster)
colors_cluster <- Polychrome::alphabet.colors(n = length(texts_cluster_uniq))
colors_cluster <- Polychrome::palette36.colors(n = length(texts_cluster_uniq))

names(colors_cluster) <- texts_cluster_uniq

## make plot
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Text_Cluster), shape = 16, alpha = 0.8, size = 1)
p <- p + scale_color_manual(values = colors_cluster, na.translate = T)
# p <- p + ggtitle(label = paste0("Tumor-cell subclusters for sample ", easyid_tmp), subtitle = "original")
p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = NULL, nrow = 2))
p <- p + theme_void()
p <- p + theme(legend.position = "bottom")
p <- p + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank())
p
## save plot
png2write <- paste0(dir_out_png, easyid_tmp, ".png")
png(filename = png2write, width = 900, height = 1000, res = 150)
print(p)
dev.off()


## make different output files
dir_out_png <- paste0(dir_out, "png", "/")
dir.create(dir_out_png)
dir_out_pdf <- paste0(dir_out, "pdf", "/")
dir.create(dir_out_pdf)
colors_all <- Polychrome::dark.colors(n = 24)

# for (easyid_tmp in "C3L-00010-T1") {
for (easyid_tmp in unique(barcode2umap_df$easy_id)) {
  scrublets_df <- barcode2scrublet_df %>%
    filter(Aliquot_WU == easyid_tmp) %>%
    filter(predicted_doublet)
  barcodes_doublet <- scrublets_df$Barcode; length(barcodes_doublet)
  
  plot_data_df <- barcode2umap_df %>%
    filter(easy_id == easyid_tmp) %>%
    filter(!(barcode_tumorcellreclustered %in% barcodes_doublet)) %>%
    mutate(Name_TumorCluster = paste0("C", id_manual_cluster_w0+1))
  cellnumber_percluster_df <- plot_data_df %>%
    select(Name_TumorCluster) %>%
    table() %>%
    as.data.frame() %>%
    rename(Name_TumorCluster = ".")
  plot_data_df <- plot_data_df %>%
    mutate(Name_TumorCluster = ifelse(Name_TumorCluster == "CNA" | Name_TumorCluster %in% cellnumber_percluster_df$Name_TumorCluster[cellnumber_percluster_df$Freq < 50], "Minor cluster (<50 cells)", Name_TumorCluster))
  
  ## make color for each cluster
  names_cluster_tmp <- sort(unique(plot_data_df$Name_TumorCluster))
  length_clusters <- length(names_cluster_tmp)
  if("Minor cluster (<50 cells)" %in% names_cluster_tmp) {
    uniq_cluster_colors <- c(colors_all[1:(length_clusters-1)], "grey40")
  } else {
    uniq_cluster_colors <-colors_all[1:length_clusters]
  }
  names(uniq_cluster_colors) <- names_cluster_tmp
  

  
  # file2write <- paste0(dir_out_pdf, easyid_tmp,".pdf")
  # pdf(file2write, width = 5, height = 5.5, useDingbats = F)
  # print(p)
  # dev.off()
}
