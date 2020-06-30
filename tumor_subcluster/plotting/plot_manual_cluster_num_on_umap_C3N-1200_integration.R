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
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200427.v1/meta_data.20200427.v1.tsv", data.table = F)
## input srat object
srat <- readRDS(file = "./Resources/snRNA_Processed_Data/Integration/C3N-01200.Tumor_Segments.Merged.20200319.v1.RDS")

# get umap info -----------------------------------------------------------
barcode2umap_df <- FetchData(object = srat, vars = c("UMAP_1", "UMAP_2", "orig.ident"))
barcode2umap_df$barcode_integration_case <- rownames(barcode2umap_df)
barcode2umap_df <- barcode2umap_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode_integration_case, pattern = "_", n = 3)[,1])

# add umap to the manual cluster id ----------------------------------------------------------------------
plot_data_df <- merge(barcode2manualsubcluster_df, barcode2umap_df,
                                     by.x = c("orig.ident", "individual_barcode"),
                                     by.y = c("orig.ident", "barcode_individual"))
plot_data_df$aliquot_wu <- mapvalues(x = plot_data_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
# plot by each aliquot ----------------------------------------------------
plot_data_df <- plot_data_df %>%
  # filter(Cell_type.shorter != "Unknown") %>%
  # filter(!is.na(Id_TumorManualCluster)) %>%
  mutate(Name_TumorCluster = paste0(aliquot_wu, "_C", (Id_TumorManualCluster+1)))
rownames(plot_data_df) <- plot_data_df$individual_barcode

## make color for each cluster
uniq_cluster_colors <- Polychrome::dark.colors(n = length(unique(plot_data_df$Name_TumorCluster)))
names(uniq_cluster_colors) <- unique(plot_data_df$Name_TumorCluster)
uniq_cluster_colors[grepl(x = names(uniq_cluster_colors), pattern = "CNA")] <- "grey80"

## make plot
p <- ggplot()
p <- p + geom_point(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = Name_TumorCluster), shape = 16, alpha = 0.8, size = 0.5)
# p <- p + ggtitle(label = paste0(aliquot_wu_tmp, " Tumor-Cell-Only Clustering"))
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
png2write <- paste0(dir_out, "C3N-01200", ".TumorCellOnAllCellIntegration.", run_id, ".png")
png(filename = png2write, width = 900, height = 1000, res = 150)
print(p)
dev.off()