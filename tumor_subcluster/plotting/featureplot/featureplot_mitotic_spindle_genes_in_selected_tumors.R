# Yige Wu @WashU Aug 2020

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
## input seurat paths
paths_srat_df <- fread(data.table = F, input = "./Data_Freezes/V1/snRNA/Tumor_Cell_Reclustered/Paths_TumorCellOnlyReclustered_SeuratObject.20201127.v1.tsv")
## input the barcode-manualsubcluster info
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20200902.v1/scrublet.run20200902_adj_cutoff.united_outputs.tsv", data.table = F)
## input genes to plot
gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/map_genes2pathway_w_avgexp_sct_data/20210414.v1/TumorCluster.DEG2Pathway.tsv")
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/unite_degs/unite_degs_for_tumor_manualcluster/20210413.v1/TumorManualCluster.DEGs.Wilcox.Minpct0.1.Logfc0.min.diff.pct0.1.tsv")

# preprocess ---------------------------------------------------
genes_filtered_df <- gene2pathway_df %>%
  filter(GeneSet_Name == "HALLMARK_MITOTIC_SPINDLE")

# process each sample -----------------------------------------------------
easyid_tmp <- "C3L-01313-T1"
easyid_tmp <- "C3L-00583-T1"

dir_out_tmp <- paste0(dir_out, easyid_tmp, "/"); dir.create(path = dir_out_tmp)
## input seurat object,
path_srat_katmai <- paths_srat_df$Path_katmai[paths_srat_df$Aliquot.snRNA.WU == easyid_tmp]
path_srat_relative <- gsub(x = path_srat_katmai, pattern = "\\/diskmnt\\/Projects\\/ccRCC_scratch\\/ccRCC_snRNA\\/", replacement = "./")
srat <- readRDS(file = path_srat_relative)
DefaultAssay(srat) <- "RNA"
## subset
scrublets_df <- barcode2scrublet_df %>%
  filter(Aliquot_WU == easyid_tmp) %>%
  filter(predicted_doublet)
barcodes_keep <- rownames(srat@meta.data)
barcodes_keep <- barcodes_keep[!(barcodes_keep %in% scrublets_df$Barcodes)]
srat <- subset(x = srat, cells = barcodes_keep)
## specify genes to plot
sampleid_tmp <- paths_srat_df$Sample[paths_srat_df$Aliquot.snRNA.WU == easyid_tmp]
degs_filtered_df <- degs_df %>%
  filter(p_val_adj < 0.001) %>%
  filter(easy_id == sampleid_tmp) %>%
  filter(gene %in% genes_filtered_df$GeneSymbol)
genes_plot <- unique(degs_filtered_df$gene)
genes_plot <- c("MKI67", genes_plot)
## write output
file2write <- paste0(dir_out_tmp, easyid_tmp, ".DEGs.tsv")
write.table(x = degs_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
# plot by gene ------------------------------------------------------------
# for (gene_plot in c("MKI67")) {
for (gene_plot in genes_plot) {
  # for (gene_plot in "AXL") {
  p <- FeaturePlot(object = srat, features = gene_plot, order = T, min.cutoff = "q10", max.cutoff = "q90", label = F)
  p <- p + scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))
  ## write output
  file2write <- paste0(dir_out_tmp, gene_plot, ".png")
  png(file2write, width = 1000, height = 1000, res = 150)
  print(p)
  dev.off()
}

