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

# get cluster 17 barcodes -------------------------------------------------
barcodes_by_aliquot_by_cluster_df <- barcode2cluster_df %>%
  rename(aliquot = orig.ident) %>%
  rename(cluster = ident) %>%
  filter(cluster == 17) %>%
  select(aliquot, cluster) %>%
  table() %>%
  as.data.frame()
barcodes_by_aliquot_by_cluster_df %>%
  arrange(-Freq)
## Looks like only CPT0075720013 have enough normal-like cells
### get the barcode in cluster 17
c17_barcodes_clean <- barcode2cluster_df %>%
  rename(aliquot = orig.ident) %>%
  rename(cluster = ident) %>%
  filter(cluster == 17) %>%
  mutate(barcode_clean = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1])

# input seurat object for CPT0075720013 -----------------------------------
seurat_obj <- readRDS("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0075720013/pf1000_fmin200_fmax10000_cmin1000_cmax10000_mito_max0.1/CPT0075720013_processed.rds")
## get dimplot data
p <- DimPlot(seurat_obj, reduction = "umap", label = T, label.size	= 5, repel = T)
label_data <- p$layers[[2]]$data
plot_data_df <- p$data 
## annotate barcodes clustered with the NAT sample (cluster 17) in the integrated dataset
plot_data_df$Is_in_normal_like_cluster <- (rownames(plot_data_df) %in% c17_barcodes_clean$barcode_clean[c17_barcodes_clean$aliquot == "CPT0075720013"])
## make new dimplot
p <- ggplot() +
  geom_point(data = plot_data_df, mapping = aes(UMAP_1, UMAP_2, color=Is_in_normal_like_cluster), alpha = 1, size = 0.3)
p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
p



