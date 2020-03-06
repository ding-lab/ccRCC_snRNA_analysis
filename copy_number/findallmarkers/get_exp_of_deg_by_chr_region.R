# Yige Wu @WashU Feb 2020
## getting the expressin of CNv-related DEGs in cells with or without CNVs per sample

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
## input the paths for individual seurat object
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input the barcode to CNV state info
barcode2cnv_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/annotate_barcode_cnv/20200228.v1/Expected_CNV_State_By_Chr_Region_By_Barcode.20200228.v1.tsv", data.table = F)
## input DEG united list
deg_by_chr_region_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/findallmarkers/unite_degs_by_chr_region/20200228.v1/FindMarkers.Wilcox.ExpectedCNV_vs_Neutral.20200228.v1.tsv", data.table = F)

# plot by chr region ------------------------------------------------------
## count the frequency of each gene associated with each chr_region
deg_count_by_chr_region_df <- deg_by_chr_region_df %>%
  select(gene, chr_region) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  arrange(desc(Freq))
## filter genes associated with the same CNV across > 1 sample
reoccurent_genes <- deg_count_by_chr_region_df %>% 
  filter(Freq > 1)

plot_data_df <- NULL
for (aliquot_tmp in unique(deg_by_chr_region_df$aliquot)) {
  ## input seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  seurat_obj_path
  srat <- readRDS(file = seurat_obj_path)
  
  ## get barcode to cnv info for this aliquot
  barcode2cnv_aliquot_df <- barcode2cnv_df %>%
    filter(aliquot == aliquot_tmp)
  
  ## set the identities to each barcodes in the meta data
  srat@meta.data$barcode <- rownames(srat@meta.data)
  metadata_tmp <- srat@meta.data
  nrow(metadata_tmp)
  metadata_tmp <- merge(metadata_tmp, barcode2cnv_aliquot_df, by = c("barcode"), all.x = T)
  nrow(metadata_tmp)
  srat@meta.data <- metadata_tmp
  
  for (chr_region_tmp in unique(deg_by_chr_region_df$chr_region[deg_by_chr_region_df$aliquot == aliquot_tmp])) {
    srat_tmp <- srat
    ## make dotplot to get the expression data
    Idents(srat_tmp) <- chr_region_tmp
    ## get genes to plot
    genes2plot_tmp <- deg_by_chr_region_df$gene[deg_by_chr_region_df$chr_region ==chr_region_tmp & deg_by_chr_region_df$aliquot == aliquot_tmp]
    genes2plot_tmp <- genes2plot_tmp[genes2plot_tmp %in% genes2plot$gene]
    genes2plot_tmp <- intersect(genes2plot_tmp, rownames(srat_tmp@assays$RNA@counts))
    if (length(genes2plot_tmp) == 0) {
      next()
    }
    p <- DotPlot(object = srat_tmp, features = genes2plot_tmp)
    plot_data_tmp <- p$data
    plot_data_tmp$aliquot <- aliquot_tmp
    plot_data_tmp$chr_region <- chr_region_tmp
    plot_data_df <- rbind(plot_data_tmp, plot_data_df)
  }
}
## save plot data
write.table(x = plot_data_df, file = paste0(dir_out, "CNV_Marker_Average_Exp_Data.", run_id, ".tsv"), sep = "\t", quote = F, row.names = F)

