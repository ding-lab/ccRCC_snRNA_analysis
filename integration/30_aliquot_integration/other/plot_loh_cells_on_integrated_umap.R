# Yige Wu @WashU Feb 2020
## for plotting the loop of Henle cells onto the UMAP

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
## input cluster to cell type table for the individual sample
cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20200207.v2.tsv", data.table = F)
## input meta data
id_metadata_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)

# gather the barcodes for samples with LOH assigned cluster ---------------
cluster2celltype_loh_df <- cluster2celltype_df %>%
  filter(Most_Enriched_Cell_Type1 == "Loop of Henle")
## loop for each aliquot
barcode_cluster_keep_df <- NULL
for (snRNA_aliquot_id_tmp in unique(cluster2celltype_loh_df$Aliquot)) {
  path_file <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/inputs/annotations_file/Individual.20200207.v1/", snRNA_aliquot_id_tmp, ".Barcode_Annotation.txt")
  barcode_cluster_tmp <- fread(input = path_file, data.table = F)
  clusters_keep <- cluster2celltype_loh_df$Cluster[cluster2celltype_loh_df$Aliquot == snRNA_aliquot_id_tmp]
  barcode_cluster_keep_tmp <- barcode_cluster_tmp %>%
    filter(V2 %in% clusters_keep) %>%
    mutate(aliquot = snRNA_aliquot_id_tmp)
  barcode_cluster_keep_df <- rbind(barcode_cluster_keep_tmp, barcode_cluster_keep_df)
}
## rename columns
barcode_cluster_keep_df <- barcode_cluster_keep_df %>%
  rename(individual_cluster = V2) %>%
  rename(individual_barcode = V1)
## reformat the barcode from integrated object
barcode2cluster_df <- barcode2cluster_df %>%
  mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]) %>%
  rename(aliquot = orig.ident)

# plot umap ---------------------------------------------------------------
plot_data_df <- merge(barcode2cluster_df, barcode_cluster_keep_df, by = c("aliquot", "individual_barcode"), all.x = T)
plot_data_df <- plot_data_df %>%
  mutate(color_cat = ifelse(is.na(individual_cluster), NA, aliquot)) %>%
  arrange(desc(color_cat))

p <- ggplot()
p <- p + geom_point(data = plot_data_df[is.na(plot_data_df$color_cat),], mapping = aes(UMAP_1, UMAP_2, color=color_cat), alpha = 1, size = 0.3)
p <- p + geom_point(data = plot_data_df[!is.na(plot_data_df$color_cat),], mapping = aes(UMAP_1, UMAP_2, color=color_cat), alpha = 1, size = 0.3)
# p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
p <- p + theme(legend.position = "top")
p
file2write <- paste(dir_out, "30_aliquot_integration", ".LOH_cells_on_umap.", run_id, ".png", sep="")
png(file2write, width = 800, height = 900, res = 150)
print(p)
dev.off()
