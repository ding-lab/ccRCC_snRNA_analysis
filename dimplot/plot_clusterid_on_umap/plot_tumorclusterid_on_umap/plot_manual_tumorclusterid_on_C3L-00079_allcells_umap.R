# Yige Wu @WashU Oct 2020

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
barcode2manualsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_C3L-00079_tumorclusterid/20201008.v1/TumorCellReclustered.BarcodeInfo.20201008.v1tsv", data.table = F)
## umap info for the all cells clustering for individual sample
umapinfo_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample/20200717.v1/Barcode2MetaData.20200717.v1.tsv")

# make data frame for plotting --------------------------------------------
plot_data_df <- umapinfo_df %>%
  filter(aliquot == "CPT0001260013") %>%
  select(individual_barcode, UMAP_1, UMAP_2)
## merge
plot_data_df <- merge(plot_data_df, 
                      barcode2manualsubcluster_df %>%
                        select(barcode_tumorcellreclustered, id_manual_cluster), 
                      by.x = c("individual_barcode"), by.y = c("barcode_tumorcellreclustered"), all.x = T)
## add cluster name
plot_data_df <- plot_data_df %>%
  mutate(ClusterName = ifelse(is.na(id_manual_cluster), "Non-tumor", paste0("C", (id_manual_cluster+1))))
table(plot_data_df$ClusterName)
# make plotting parameters ------------------------------------------------
## make distinguishable colors
colors_cluster <- RColorBrewer::brewer.pal(n = length(table(plot_data_df$ClusterName)), name = "Dark2")
names(colors_cluster) <- names(table(plot_data_df$ClusterName))

# Dimplot -----------------------------------------------------------------
## plot
## make plot
p <- ggplot()
p <- p + geom_point(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = ClusterName), shape = 16, alpha = 0.8, size = 0.5)
p <- p + scale_color_manual(values = colors_cluster, na.translate = T)
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

# write output ------------------------------------------------------------
file2write <- paste0(dir_out,  "C3L-00079", ".png")
png(file2write, width = 900, height = 800, res = 150)
print(p)
dev.off()
