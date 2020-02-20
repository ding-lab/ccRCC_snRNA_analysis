# Yige Wu @WashU Feb 2020
## running on local
## for plotting the fraction of cells belong to different immune cell types in each aliquot in the integrated immune dataset

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
barcode2cluster_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_integrated_data/fetch_data/20200214.v1/integration.202002012.v3.immune_reclustered.20200213.v2.umap_data.20200214.v1.tsv", data.table = F)
## input cluster to cell type table
cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/integration.immune.cluster2celltype.20200214.v1.csv", data.table = F)
## input meta data
id_metadata_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)

# transform and merge data --------------------------------------------------------------
## for each aliquot count the total number of barcodes
barcodes_by_aliquot_df <- barcode2cluster_df %>%
  select(orig.ident) %>%
  table() %>%
  as.data.frame() %>%
  rename(aliquot = '.')
## for each aliquot, count the barcodes for each cluster
barcodes_by_aliquot_by_cluster_df <- barcode2cluster_df %>%
  select(orig.ident, ident) %>%
  rename(aliquot = orig.ident) %>%
  rename(cluster = ident) %>%
  table() %>%
  as.data.frame()

## merge barcode counts
barcodes_by_aliquot_by_cluster_df <- base::merge(barcodes_by_aliquot_by_cluster_df, barcodes_by_aliquot_df, 
                                                 by = c("aliquot"),
                                                 suffixes = c("_by_aliquot_by_cluster", "_by_aliquot"), all.x = T)
## get fraction
barcodes_by_aliquot_by_cluster_df <- barcodes_by_aliquot_by_cluster_df %>%
  mutate(Frac_by_aliquot_by_cluster = Freq_by_aliquot_by_cluster/Freq_by_aliquot)

## merge cluster2celltype
barcodes_by_aliquot_by_cluster_df <- merge(barcodes_by_aliquot_by_cluster_df, cluster2celltype_df,
                                           by.x = c("cluster"), by.y = c("Cluster"),
                                           all.x = T)

## map aliquot to case
barcodes_by_aliquot_by_cluster_df$case <- mapvalues(x = barcodes_by_aliquot_by_cluster_df$aliquot, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

## make a column for cell type to show
celltypes2show <- as.vector(barcodes_by_aliquot_by_cluster_df$Most_Enriched_Cell_Type4)
celltypes2show[celltypes2show == ""] <- as.vector(barcodes_by_aliquot_by_cluster_df$Most_Enriched_Cell_Type3[celltypes2show == ""])
celltypes2show[celltypes2show == ""] <- as.vector(barcodes_by_aliquot_by_cluster_df$Most_Enriched_Cell_Type2[celltypes2show == ""])
celltypes2show[celltypes2show == ""] <- as.vector(barcodes_by_aliquot_by_cluster_df$Most_Enriched_Cell_Type1[celltypes2show == ""])
celltypes2show[celltypes2show == ""] <- as.vector(barcodes_by_aliquot_by_cluster_df$Most_Enriched_Cell_Group[celltypes2show == ""])
barcodes_by_aliquot_by_cluster_df$Cell_Type2Show <- celltypes2show

# make barplot by cell type fraction------------------------------------------------------------
plot_df <- barcodes_by_aliquot_by_cluster_df
plot_df <- plot_df %>%
  arrange(desc(aliquot))
aliquots <- unique(plot_df)
plot_df$aliquot <- factor()
## make color palette
unique(barcodes_by_aliquot_by_cluster_df$Cell_Type2Show)
celltype2colors <- c("Macrophages M1" = "#a6cee3",
                     "CD4+ T-cells" = "#1f78b4",
                     "CD8+ T-cells" = "#b2df8a",
                     "Macrophages M2" = "#33a02c",
                     "Unknown" = "grey50",
                     "CD8+ T-cells & CD4+ T-cells" = "#fb9a99",
                     "cDC1" = "#e31a1c",
                     "B-cells" = "#fdbf6f",
                     "Plasma cells" = "#ff7f00",
                     "Myeloid lineage immune cells" = "#cab2d6",
                     "Macrophages M1&M2" = "#6a3d9a",
                     "NK-cells" = "#ffff99")

#b15928
p <- ggplot()
p <- p + geom_col(data = plot_df, mapping = aes(x = aliquot, y = Frac_by_aliquot_by_cluster, fill = Cell_Type2Show), color = "black", position = "stack")
p <- p + scale_fill_manual(values = celltype2colors)
p <- p + facet_grid(.~case, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5),
               strip.text.x = element_text(angle = 90, vjust = 0.5, size = 12, face = "bold"),
               panel.border = element_rect(color = "black", fill = NA, size = 1),
               strip.placement = "outside")
## save barplot
file2write <- paste0(dir_out, "integration.202002012.v3.immune_reclustered.20200213.v2.Immune_CellType_Fraction.by_aliquot_by_cluster.", run_id, ".png")
png(file = file2write, width = 3000, height = 1500, res = 150)
print(p)
dev.off()


