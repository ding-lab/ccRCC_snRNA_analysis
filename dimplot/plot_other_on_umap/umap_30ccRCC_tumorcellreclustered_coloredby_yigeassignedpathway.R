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
integrated_umap_df <- fread(input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/annotate_clusterid_on_30ccRCCtumorcellreclustered_byresolution/20220406.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.20220406.v1.tsv", data.table = F)
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210805.v1/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

# make plot data----------------------------------------------------------
integrated_umap_df$easy_id <- mapvalues(x = integrated_umap_df$orig.ident, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Aliquot.snRNA.WU))
integrated_umap_df <- integrated_umap_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1])
plot_data_df <- merge(x = integrated_umap_df, y = barcode2tumorsubcluster_df, 
                         by.x = c("easy_id", "barcode_individual"),
                         by.y = c("easy_id", "barcode"),
                         all.x = T)
enrich_df <- enrich_df %>%
  mutate(cluster_name.merge = gsub(x = cluster_name, pattern = "\\.", replacement = "-"))
plot_data_df <- merge(x = plot_data_df, y = enrich_df %>%
                        select(cluster_name.merge, Cell_cycle, EMT, Immune, mTOR), 
                      by.x = c("Cluster_Name"),
                      by.y = c("cluster_name.merge"),
                      all.x = T)

# make plots --------------------------------------------------------------
genemodule_tmp <- "mTOR"
for (genemodule_tmp in c("Cell_cycle", "Immune", "EMT", "mTOR")) {
  plot_data_df[, "cell_group"] <- ifelse(is.na(plot_data_df[, genemodule_tmp]), "Not assigned",
                                         ifelse(plot_data_df[, genemodule_tmp], paste0(genemodule_tmp, " enriched (original)"), "Not enriched"))
  plot_data_df <- plot_data_df %>%
    arrange(factor(cell_group, levels = rev(c(paste0(genemodule_tmp, " enriched (original)"),
                                          "Not enriched",
                                          "Not assigned"))))
  colors_cellgroup <- c("red", "grey80", "grey50")
  names(colors_cellgroup) <- c(paste0(genemodule_tmp, " enriched (original)"),
                               "Not enriched",
                               "Not assigned")
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df, 
                           mapping = aes(x = UMAP_1, y = UMAP_2, color = cell_group),
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
  p <- p + theme(legend.position="bottom", aspect.ratio=1)
  p
  # ## save as pdf
  # file2write <- paste0(dir_out, "cellgroup_on_umap.", ".pdf")
  # pdf(file = file2write, width = 8, height = 9, useDingbats = F)
  # print(p)
  # dev.off()
  ## save as png
  file2write <- paste0(dir_out, genemodule_tmp, "_on_umap.", ".png")
  png(filename = file2write, width = 800, height = 850, res = 150)
  print(p)
  dev.off()
}
