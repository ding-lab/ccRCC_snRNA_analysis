# Yige Wu @WashU Sep 2021

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
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_individual_tumorcellreclustered_on_katmai/20210805.v1/MetaData_TumorCellOnlyReclustered.20210805.v1.tsv")
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)
epi_group_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_epithelial_group_bytumorcluster/20211011.v1/Tumorcluster_EpithelialGroup.20211011.v1.tsv")

# merge data --------------------------------------------------------------
barcode_merged_df <- merge(x = barcode2umap_df, y = barcode2tumorsubcluster_df, 
                         by.x = c("easy_id", "barcode_tumorcellreclustered"),
                         by.y = c("easy_id", "barcode"),
                         all.x = T)
## add enrichment type
cluster2group_df <- epi_group_df %>%
  mutate(sample = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
  mutate(sample = gsub(x = sample, pattern = "\\.", replacement = "-")) %>%
  mutate(cluster = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,2]) %>%
  mutate(cluster_name2 = paste0(sample, "_", cluster))

table(barcode_merged_df$Cluster_Name)
barcode_merged_df <- merge(x = barcode_merged_df, y = cluster2group_df, by.x = c("Cluster_Name"), by.y = c("cluster_name2"), all.x = T)

# plot by each aliquot ----------------------------------------------------
## make different output files
dir_out_png <- paste0(dir_out, "png", "/")
dir.create(dir_out_png)
dir_out_pdf <- paste0(dir_out, "pdf", "/")
dir.create(dir_out_pdf)
# colors_enrich_type <- RColorBrewer::brewer.pal(n = 8, name = "Set2")[c(1,6,3,4,8)]
## make colors for cell type
colors_tumorgroup_sim <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(1, 5, 4, 2)],  colors_cellgroup14[c("Normal epithelial cells")])
names(colors_tumorgroup_sim) <- c("EMT",  "Epi-L", "Epi-M", "Epi-H", "PT")

# for (aliquot_tmp in "C3L-00416-T1") {
for (aliquot_tmp in unique(barcode_merged_df$easy_id)) {
  scrublets_df <- barcode2scrublet_df %>%
    filter(Aliquot_WU == aliquot_tmp) %>%
    filter(predicted_doublet)
  
  plot_data_df <- barcode_merged_df %>%
    filter(easy_id == aliquot_tmp) %>%
    filter(!(barcode_tumorcellreclustered %in% scrublets_df$Barcodes))
  
  ## make color for each cluster

  ## make png
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = epithelial_group), shape = 16, alpha = 1, size = 1)
  p <- p + scale_color_manual(values = colors_tumorgroup_sim, na.translate = T)
  p <- p + theme_bw()
  p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
  p <- p + theme_void()
  p <- p + theme(legend.position = "none")
  p <- p + ggtitle(label = str_split_fixed(string = plot_data_df$cluster_name.figure[1], pattern = "_", n = 2)[,1])
  p <- p + theme(title = element_text(size = 17))
  p <- p + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank())
  # save plot
  png2write <- paste0(dir_out_png, aliquot_tmp, ".png")
  png(filename = png2write, width = 900, height = 900, res = 150)
  print(p)
  dev.off()
  
  # p <- ggplot()
  # p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = epithelial_group), shape = 16, alpha = 1, size = 1)
  # p <- p + scale_color_manual(values = colors_tumorgroup_sim, na.translate = T)
  # p <- p + theme_bw()
  # p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
  # p <- p + theme_void()
  # p <- p + theme(legend.position = "none")
  # p <-  p + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))
  # p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #                panel.background = element_blank())
  # file2write <- paste0(dir_out_pdf, aliquot_tmp,".pdf")
  # pdf(file2write, width = 5, height = 5, useDingbats = F)
  # print(p)
  # dev.off()
}
