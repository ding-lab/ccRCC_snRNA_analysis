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
## set genes to plot
genes2plot <- c("VHL", "EPAS1", "HIF1A", "PFKFB4")

# input srat object and edit meta data ------------------------------------
## add case id
barcode2manualsubcluster_df$id_case <- mapvalues(x = barcode2manualsubcluster_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
barcode2manualsubcluster_df$id_aliquot_wu <- mapvalues(x = barcode2manualsubcluster_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## get barcode to cnv info for this aliquot
barcode2manualsubcluster_aliquot_df <- barcode2manualsubcluster_df %>%
  filter(id_case == "C3N-01200") %>%
  filter(!is.na(Id_TumorManualCluster)) %>%
  mutate(Name_Cluster = paste0(id_aliquot_wu, "_C", (Id_TumorManualCluster+1)))
unique(barcode2manualsubcluster_aliquot_df$Name_Cluster)

# get umap info -----------------------------------------------------------
metadata_df <- srat@meta.data
metadata_df$barcode_integrated_case <- rownames(metadata_df)
metadata_df <- metadata_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode_integrated_case, pattern = "_", n = 3)[,1])
barcode2manualsubcluster_aliquot_df <- merge(barcode2manualsubcluster_aliquot_df, metadata_df, 
                                             by.x = c("orig.ident", "individual_barcode"), by.y = c("orig.ident", "barcode_individual"), all.x = T)

## filter down to only tumor cells
srat_plot <- subset(x = srat, cells = barcode2manualsubcluster_aliquot_df$barcode_integrated_case)
## change meta data
metadata_df <- srat_plot@meta.data
metadata_df$Name_Cluster <- mapvalues(x = rownames(metadata_df), from = barcode2manualsubcluster_aliquot_df$barcode_integrated_case, to = as.vector(barcode2manualsubcluster_aliquot_df$Name_Cluster))
srat_plot@meta.data <- metadata_df
### set the identities to cluster in the meta data
Idents(object = srat_plot) <- "Name_Cluster"


# plot --------------------------------------------------------------------
for (gene_tmp in genes2plot) {
  p <- VlnPlot(object = srat_plot, features = gene_tmp, group.by = "Name_Cluster", pt.size = 0)
  p <- p + theme(axis.text.x = element_text(angle = 90))
  p <- p + theme(legend.position = "none")
  png(filename = paste0(dir_out, "C3N-01200", ".", gene_tmp, ".png"),width = 700, height = 500, res = 150)
  print(p)
  dev.off()
}
