# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(monocle)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
## input paths to the monocle objects
obj_monocle <- readRDS(file = "./Resources/snRNA_Processed_Data/Monocle/outputs/CellTypeVer.20200828.v1/C3L-00079_PooledPT_SelfFib_ByCellType/combined_subset_pseudotime_qval_1e-10.rds")
## input meta data
barcode2manualsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_C3L-00079_tumorclusterid/20201008.v1/TumorCellReclustered.BarcodeInfo.20201008.v1tsv", data.table = F)

# plot ---------------------------------------------------------
pdata_df <- as.data.frame(pData(obj_monocle))
pdata_df <- pdata_df %>%
  mutate(id_cell = paste0(orig.ident, "_", original_barcode))
barcode2manualsubcluster_df <- barcode2manualsubcluster_df %>%
  mutate(id_cell = paste0(orig.ident, "_", barcode_tumorcellreclustered)) %>%
  mutate(Name_Cluster = paste0("C", (id_manual_cluster+1)))
## map tumor subcluster id
pdata_df$Name_Cluster <- mapvalues(x = pdata_df$id_cell, from = barcode2manualsubcluster_df$id_cell, to = as.vector(barcode2manualsubcluster_df$Name_Cluster))
pdata_df$Name_Cluster[pdata_df$Name_Cluster == pdata_df$id_cell] <- "Non-tumor"
table(pdata_df$Name_Cluster)
pData(obj_monocle)$Name_Cluster <- pdata_df$Name_Cluster
## specify cluster number
colors_clustername <- c(RColorBrewer::brewer.pal(n = 5, name = "Dark2"), "grey50")
names(colors_clustername) <- names(table(pdata_df$Name_Cluster))
swatch(colors_clustername)

## plot
p <- plot_cell_trajectory(obj_monocle, color_by = "Name_Cluster",cell_size=0.6)
p <- p + scale_color_manual(values = colors_clustername)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(aspect.ratio=1)
## write output
file2write <- paste0(dir_out, "seuratclusters_on_trajecctory", ".png")
png(file2write, width = 1000, height = 1100, res = 150)
print(p)
dev.off()
