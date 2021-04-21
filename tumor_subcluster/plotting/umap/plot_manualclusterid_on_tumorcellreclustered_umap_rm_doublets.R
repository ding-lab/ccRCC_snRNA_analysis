# Yige Wu @WashU Dec 2020

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
## input barcode2UMAP
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_individual_tumorcellreclustered_on_katmai/20201127.v1/MetaData_TumorCellOnlyReclustered.20201127.v1.tsv")
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20201130.v1/Barcode2TumorSubclusterId.20201130.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20200902.v1/scrublet.run20200902_adj_cutoff.united_outputs.tsv", data.table = F)

# merge data --------------------------------------------------------------
barcode2umap_df <- merge(x = barcode2umap_df, y = barcode2tumorsubcluster_df, 
                         by.x = c("easy_id", "barcode_tumorcellreclustered"),
                         by.y = c("easy_id", "barcode"),
                         all.x = T)


# plot by each aliquot ----------------------------------------------------
## make different output files
dir_out_png <- paste0(dir_out, "png", "/")
dir.create(dir_out_png)
dir_out_pdf <- paste0(dir_out, "pdf", "/")
dir.create(dir_out_pdf)

for (aliquot_tmp in "C3L-00416-T1") {
# for (aliquot_tmp in unique(barcode2umap_df$easy_id)) {
  scrublets_df <- barcode2scrublet_df %>%
    filter(Aliquot_WU == aliquot_tmp) %>%
    filter(predicted_doublet)
  
  plot_data_df <- barcode2umap_df %>%
    filter(easy_id == aliquot_tmp) %>%
    filter(!(barcode_tumorcellreclustered %in% scrublets_df$Barcodes)) %>%
    mutate(Name_TumorCluster = paste0("C", id_manual_cluster_w0+1))

  
  ## make color for each cluster
  uniq_cluster_colors <- Polychrome::dark.colors(n = length(unique(plot_data_df$Name_TumorCluster)))
  names(uniq_cluster_colors) <- unique(plot_data_df$Name_TumorCluster)
  
  ## make plot
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Name_TumorCluster), shape = 16, alpha = 0.8, size = 1)
  p <- p + scale_color_manual(values = uniq_cluster_colors, na.translate = T)
  p <- p + theme_bw()
  p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
  p <- p + theme_void()
  p <- p + theme(legend.position = "top")
  p <-  p + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank())
  ## save plot
  png2write <- paste0(dir_out_png, aliquot_tmp, ".png")
  png(filename = png2write, width = 900, height = 1000, res = 150)
  print(p)
  dev.off()
  
  # file2write <- paste0(dir_out_pdf, aliquot_tmp,".pdf")
  # pdf(file2write, width = 5, height = 5.5, useDingbats = F)
  # print(p)
  # dev.off()
}
