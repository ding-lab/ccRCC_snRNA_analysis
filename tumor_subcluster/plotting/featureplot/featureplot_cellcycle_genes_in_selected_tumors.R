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
## input genes to plot
gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/map_genes2pathway_w_avgexp_sct_data/20210414.v1/TumorCluster.DEG2Pathway.tsv")
## input seurat paths
paths_srat_df <- fread(data.table = F, input = "./Data_Freezes/V1/snRNA/Tumor_Cell_Reclustered/Paths_TumorCellOnlyReclustered_SeuratObject.20201127.v1.tsv")
## input the barcode-manualsubcluster info
barcode2subclusterid_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20201130.v1/Barcode2TumorSubclusterId.20201130.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20200902.v1/scrublet.run20200902_adj_cutoff.united_outputs.tsv", data.table = F)
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210421.v3/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

# preprocess ---------------------------------------------------
## specify genes to plot
genes_plot <- c("ABR", "MYO9B", "NIN", "CNTRL", "CAPZB")

## get samples to plot
enrich_filtered_df <- enrich_df %>%
  filter(Cell_cycle) %>%
  mutate(easyid = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
  mutate(easyid = gsub(x = easyid, pattern = "\\.", replacement = "-"))

# process by sample -------------------------------------------------------
easyid_tmp <- "C3L-00416-T1"
easyid_tmp <- "C3L-00096-T1"
for (easyid_tmp in unique(enrich_filtered_df$easyid)) {
  sampleid_tmp <- paths_srat_df$Sample[paths_srat_df$Aliquot.snRNA.WU == easyid_tmp]
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
  
  # plot by gene ------------------------------------------------------------
  for (gene_plot in genes_plot) {
    dir_out_tmp <- paste0(dir_out, gene_plot, "/"); dir.create(dir_out_tmp)
    file2write <- paste0(dir_out_tmp, easyid_tmp, "_", gene_plot, ".png")
    if (!file.exists(file2write)) {
      # for (gene_plot in "AXL") {
      p <- FeaturePlot(object = srat, features = gene_plot, order = T, min.cutoff = "q10", max.cutoff = "q90", label = F)
      p <- p + scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))
      ## write output
      png(file2write, width = 1000, height = 1000, res = 150)
      print(p)
      dev.off()
    }
  }
}
