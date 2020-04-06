# Yige Wu @WashU March 2020
## run DEG analysis for cells with certain CNV vs CN neutral

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
## set chr region to examine
chr_region_tmp <- "5q"
## input the paths for individual seurat object
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input the barcode to CNV state info
barcode2cnv_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/annotate_barcode_cnv/20200228.v1/Expected_CNV_State_By_Chr_Region_By_Barcode.20200228.v1.tsv", data.table = F)
## input 5q DEGs
deg_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/findallmarkers/examine_5q_degs/20200302.v1/5q.FindMarkers.Wilcox.ExpectedCNV_vs_Neutral.20200302.v1.tsv", data.table = F)
## input gene to pathway annotatin
gene2pathway_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/findallmarkers/examine_5q_degs/20200302.v1/5q.FindMarkers.Wilcox.ExpectedCNV_vs_Neutral.Gene2Pathway.20200302.v1.tsv", data.table = F)
## input case CNV profile
bulk_bicseq_cnv_state_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/bulk/write_sample_bicseq_cnv_profile/20200227.v1/Bulk_WGS_Chr_CNV_Profile.20200227.v1.tsv", data.table = F)
bulk_bicseq_cnv_state_df$`5q`[bulk_bicseq_cnv_state_df$Case == "C3N-00242"] <- "gain"
## filter cases by bulk WGS CNV profile
cases_expected_cnv <- bulk_bicseq_cnv_state_df[bulk_bicseq_cnv_state_df$`5q` == "gain",]
## input meta data
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
# srat_paths$Case <- mapvalues(x = srat_paths$Aliquot, from = meta_tab$Aliquot.snRNA, to = as.vector(meta_tab$Case))

# test --------------------------------------------------------------------
aliquot_tmp <- "CPT0020120013"
pathway_tmp <- "NRF2 pathway"
genes_in_pathway_tmp <- gene2pathway_df$gene[gene2pathway_df$pathway == pathway_tmp]
genes_in_pathway_tmp
deg_df.filtered <- deg_df %>%
  filter(gene %in% genes_in_pathway_tmp)

for (aliquot_tmp in unique(deg_df$aliquot))  {
  ## input seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  seurat_obj_path
  srat <- readRDS(file = seurat_obj_path)
  
  ## get barcode to cnv info for this aliquot
  barcode2cnv_aliquot_df <- barcode2cnv_df %>%
    filter(aliquot == aliquot_tmp)
  
  ### set the identities to each barcodes in the meta data
  srat@meta.data$barcode <- rownames(srat@meta.data)
  metadata_tmp <- srat@meta.data
  nrow(metadata_tmp)
  metadata_tmp <- merge(metadata_tmp, barcode2cnv_aliquot_df, by = c("barcode"), all.x = T)
  nrow(metadata_tmp)
  rownames(metadata_tmp) <- metadata_tmp$barcode
  srat@meta.data <- metadata_tmp
  
  genes2plot <- intersect(rownames(srat@assays$RNA@counts), as.vector(deg_df$gene))
  # genes2plot <- intersect(rownames(srat@assays$RNA@counts), as.vector(gene2pathway_df$gene[gene2pathway_df$pathway == pathway_tmp]))
  genes2plot
  
  Idents(object = srat) <- chr_region_tmp
  p <- DoHeatmap(object = srat, features = genes_in_pathway_tmp, raster = T)
  png(filename = paste0(dir_out, aliquot_tmp, ".", pathway_tmp, ".", run_id, ".png"),width = 1200, height = 300, res = 150)
  print(p)
  dev.off()
  
  stop("")
}





