# Yige Wu @WashU March 2020
## for each individual sample tumor cell reclustered, plot UMAP

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
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
## input barcode-to-manual grouping 
barcode2manualsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_barcode_with_manual_tumorsubcluster_id/20200616.v1/Barcode2TumorSubclusterId.20200616.v1.tsv", data.table = F)
## input barcode2UMAP
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_barcode_with_seurat_tumorsubcluster_id/20200618.v1/Barcode2SeuratClusterID.20200618.v1.tsv")
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200427.v1/meta_data.20200427.v1.tsv", data.table = F)

# plot by each aliquot ----------------------------------------------------
## add umap to the manual cluster id
barcode2manualsubcluster_df <- merge(barcode2manualsubcluster_df, barcode2umap_df,
                                     by.x = c("orig.ident", "individual_barcode"),
                                     by.y = c("orig.ident", "barcode"), all.x = T)
for (aliquot_tmp in unique(barcode2manualsubcluster_df$orig.ident)) {
  plot_data_df <- barcode2manualsubcluster_df %>%
    filter(orig.ident == aliquot_tmp) %>%
    filter(Cell_type.shorter != "Unknown") %>%
    filter(!is.na(Id_TumorManualCluster)) %>%
    mutate(Name_TumorCluster = paste0("C", (Id_TumorManualCluster+1)))
  rownames(plot_data_df) <- plot_data_df$individual_barcode
  
  ## make color for each cluster
  uniq_cluster_colors <- Polychrome::dark.colors(n = length(unique(plot_data_df$Name_TumorCluster)))
  names(uniq_cluster_colors) <- unique(plot_data_df$Name_TumorCluster)
  
  ## get aliquot id to show
  aliquot_wu_tmp <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquot_tmp]
  
  ## make plot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Name_TumorCluster), shape = 16, alpha = 0.8, size = 0.5)
  p <- p + ggtitle(label = paste0(aliquot_wu_tmp, " Tumor-Cell-Only Clustering"))
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
  png2write <- paste0(dir_out, aliquot_wu_tmp, ".TumorCellOnlyClustering.", run_id, ".png")
  png(filename = png2write, width = 900, height = 1000, res = 150)
  print(p)
  dev.off()
}
