# Yige Wu @WashU Feb 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the paths for reclustered seurat objects
srat_paths <- fread(input = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input the barcode-manualsubcluster info
barcode2subclusterid_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_barcode_with_manual_tumorsubcluster_id/20200616.v1/Barcode2TumorSubclusterId.20200616.v1.tsv", data.table = F)
## input the tumor subcluster pairs 
cnvfraction_pair_filtered_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/summarize_tumormanualsubcluster_pairs_w_cnv_diff/20200622.v1/CNVFractionDifference_between_ManualSubclusterPairs_PerGene_PerTumor.Filtered20200622.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")

# set parameters for findmarkers ------------------------------------------
## set min.pct
# min.pct.wilcox <- 0.25
min.pct.wilcox <- 0.1
logfc.threshold.wilcox <- 0.1
min.pct.wilcox <- 0
logfc.threshold.wilcox <- 0

# add additional info -------------------------
## format the barcode table
barcode2subclusterid_df$id_aliquot_wu <- mapvalues(x = barcode2subclusterid_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
barcode2subclusterid_df <- barcode2subclusterid_df %>%
  mutate(id_aliquot_cluster = paste0(id_aliquot_wu, "_C", (Id_TumorManualCluster + 1)))
barcode2subclusterid_df$id_aliquot_cluster %>% unique()
# for each aliquot, input seurat object and fetch data and write data --------------------
markers_wilcox_df <- NULL
for (aliquot_wu_tmp in unique(cnvfraction_pair_filtered_df$aliquot.wu)) {
  aliquot_tmp <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU == aliquot_wu_tmp]
  
  # input seurat object -----------------------------------------------------
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  seurat_obj_path
  srat <- readRDS(file = seurat_obj_path)
  
  # change metadata ---------------------------------------------------------
  ## change ident to manual subcluster
  ### get the barcode-subclusterid for this aliquot
  metadata_new_df <- barcode2subclusterid_df %>%
    filter(orig.ident == aliquot_tmp)
  rownames(metadata_new_df) <- metadata_new_df$barcode
  ### change meta data
  srat@meta.data <- metadata_new_df
  ### change ident
  Idents(srat) <- "id_aliquot_cluster"
  
  # get which clusters to compare -------------------------------------------
  aliquot_pair_filtered_df <- cnvfraction_pair_filtered_df %>%
    filter(aliquot.wu == aliquot_wu_tmp) %>%
    select(id_aliquot_cluster.1, id_aliquot_cluster.2) %>%
    unique()
  
  for (i in 1:nrow(aliquot_pair_filtered_df)) {
    id_aliquot_cluster.1 <- aliquot_pair_filtered_df[i, "id_aliquot_cluster.1"]
    id_aliquot_cluster.2 <- aliquot_pair_filtered_df[i, "id_aliquot_cluster.2"]
    
    ## find all markers using Wilcox testing
    markers_wilcox_tmp <- FindMarkers(object = srat, test.use = "wilcox", only.pos = F, min.pct = min.pct.wilcox, logfc.threshold = logfc.threshold.wilcox, return.thresh = 0.05, 
                                      ident.1 = id_aliquot_cluster.1, ident.2 = id_aliquot_cluster.2, verbose = T)
    
    ## add the current DEGs into the super table
    if (nrow(markers_wilcox_tmp) > 0) {
      markers_wilcox_tmp$deg_gene_symbol <- rownames(markers_wilcox_tmp)
      markers_wilcox_tmp$id_aliquot <- aliquot_tmp
      markers_wilcox_tmp$id_aliquot_wu <- aliquot_wu_tmp
      markers_wilcox_tmp$ident.1 <- id_aliquot_cluster.1
      markers_wilcox_tmp$ident.2 <- id_aliquot_cluster.2
      markers_wilcox_df <- rbind(markers_wilcox_df, markers_wilcox_tmp)
    }
  }
  stop("test")
}

# write output ------------------------------------------------------------
write.table(markers_wilcox_df, file = paste0(dir_out, "Tumormanualsubcluster.FindMarkers.Cluster_CNVLow_vs_CNVHigh.Wilcox.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox,".tsv"), 
            quote = F, sep = "\t", row.names = F)

# write.table(markers_wilcox_df%>%
#             filter(cna_gene_symbol == deg_gene_symbol), file = paste0(dir_out, "Tumormanualsubcluster.FindMarkers.Cluster_CNVLow_vs_CNVHigh.Wilcox.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox,".CNV_Genes.tsv"), 
#             quote = F, sep = "\t", row.names = F)

