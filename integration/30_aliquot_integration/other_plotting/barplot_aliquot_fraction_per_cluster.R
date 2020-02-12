# Yige Wu @WashU Feb 2020
## running on local
## for plotting the fraction of cells belong to different aliquot in each cluster in the integrated dataset

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode-aliquot-cluster table
barcode2cluster_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/20200212.v3/30_aliquot_integration.20200212.v3.umap_data.tsv", data.table = F)
## input cluster to cell type table
cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/30_aliquot_integration.allcluster2celltype.20200212.v1.tsv", data.table = F)
## input meta data
id_metadata_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)

# transform and merge data --------------------------------------------------------------
## for each cluster count the total number of barcodes
barcodes_by_cluster_df <- barcode2cluster_df %>%
  select(ident) %>%
  table() %>%
  as.data.frame() %>%
  rename(cluster = '.')
## by each cluster, count the barcodes for each aliquot
barcodes_by_aliquot_by_cluster_df <- barcode2cluster_df %>%
  select(orig.ident, ident) %>%
  rename(aliquot = orig.ident) %>%
  rename(cluster = ident) %>%
  table() %>%
  as.data.frame()
## merge barcode counts
barcodes_by_aliquot_by_cluster_df <- base::merge(barcodes_by_aliquot_by_cluster_df, barcodes_by_cluster_df, 
                                                 by = c("cluster"),
                                                 suffixes = c("_by_aliquot_by_cluster", "_by_cluster"), all.x = T)
## get fraction
barcodes_by_aliquot_by_cluster_df <- barcodes_by_aliquot_by_cluster_df %>%
  mutate(Frac_by_aliquot_by_cluster = Freq_by_aliquot_by_cluster/Freq_by_cluster)

## merge cluster2celltype
barcodes_by_aliquot_by_cluster_df <- merge(barcodes_by_aliquot_by_cluster_df, cluster2celltype_df,
                                           by.x = c("cluster"), by.y = c("Cluster"),
                                           all.x = T)

## map aliquot to case
barcodes_by_aliquot_by_cluster_df$case <- mapvalues(x = barcodes_by_aliquot_by_cluster_df$aliquot, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

# make barplot by all aliquots------------------------------------------------------------
plot_df <- barcodes_by_aliquot_by_cluster_df
p <- ggplot()
p <- p + geom_col(data = plot_df, mapping = aes(x = cluster, y = Frac_by_aliquot_by_cluster, fill = aliquot), color = "black", position = "stack")
p <- p + facet_grid(.~Most_Enriched_Cell_Group, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5),
               strip.text.x = element_text(angle = 0, vjust = 0.5, size = 12, face = "bold"),
               panel.border = element_rect(color = "black", fill = NA, size = 1),
               strip.placement = "outside")

## save barplot
file2write <- paste0(dir_out, "Cell_Group_Composition.by_aliquot_by_cluster.30_aliquot_integration.20200211.v3.", ".png")
png(file = file2write, width = 3000, height = 1500, res = 150)
print(p)
dev.off()

# make barplot by case------------------------------------------------------------
case_tmp <- "C3L-00088"
for (case_tmp in c("C3L-00088", "C3N-01200", "C3N-00733", "C3L-00416")) {
  plot_df <- barcodes_by_aliquot_by_cluster_df
  ## make aliquots from other cases the same
  tmp <- plot_df$aliquot
  tmp[plot_df$case != case_tmp] <- "other"
  plot_df$aliquot <- tmp
  ## make ggplot
  p <- ggplot()
  p <- p + geom_col(data = plot_df, mapping = aes(x = cluster, y = Frac_by_aliquot_by_cluster, fill = aliquot), color = "black", position = "stack")
  p <- p + facet_grid(.~Most_Enriched_Cell_Group, scales = "free", space = "free", drop = T)
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5),
                 strip.text.x = element_text(angle = 0, vjust = 0.5, size = 12, face = "bold"),
                 panel.border = element_rect(color = "black", fill = NA, size = 1),
                 strip.placement = "outside")
  file2write <- paste0(dir_out, "Cell_Group_Composition.by_aliquot_by_cluster.30_aliquot_integration.20200211.v3.", case_tmp,"_Aliquots.png")
  png(file = file2write, width = 3000, height = 1500, res = 150)
  print(p)
  dev.off()
}
