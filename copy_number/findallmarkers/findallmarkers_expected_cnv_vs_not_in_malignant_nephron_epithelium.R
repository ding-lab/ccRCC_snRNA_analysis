# Yige Wu @WashU Feb 2020
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
## input the paths for individual seurat object
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input the barcode to CNV state info
barcode2cnv_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/annotate_barcode_cnv/20200228.v1/Expected_CNV_State_By_Chr_Region_By_Barcode.20200228.v1.tsv", data.table = F)
## chromosome regions from which the CNVs are annotated here
chr_regions2process <- unique(ccrcc_cna_genes_df$chr_region)
chr_regions2process <- as.vector(chr_regions2process)

# for each aliquot, input seurat object and fetch data and write data --------------------
## initiate the vector to record the file path
aliquot_vec <- NULL
chr_region_vec <- NULL
deg_filepath_vec <- NULL
for (aliquot_tmp in unique(srat_paths$Aliquot)) {
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
  srat@meta.data <- metadata_tmp
  
  ## make output directory for this aliquot
  dir_out_sub1 <- paste0(dir_out, aliquot_tmp, "/")
  dir.create(dir_out_sub1)
  
  for (chr_region_tmp in chr_regions2process) {
    if (all(is.na(barcode2cnv_aliquot_df[,chr_region_tmp]))) {
      next()
    }
    srat_tmp <- srat
    Idents(srat_tmp) <- chr_region_tmp
    
    ### find all markers using Wilcox testing
    markers_wilcox <- tryCatch(expr = FindMarkers(object = srat_tmp, ident.1 = "Expected", ident.2 = "Neutral", logfc.threshold = 0.1,
                                                  only.pos = FALSE, test.use = "wilcox"),
                               error = function(e) {warning("Marker identification failed.");return(NULL)})
    if (is.null(markers_wilcox)) {
      next()
    }
    markers_wilcox$gene <- rownames(markers_wilcox)
    ### write table
    file2write <- paste0(dir_out_sub1, aliquot_tmp, ".", chr_region_tmp, ".", "ExpectedCNV_vs_Neutral", ".", ".FindAllMarkers.Wilcox.Pos", ".", run_id, ".tsv")
    write.table(markers_wilcox, file = file2write, quote = F, sep = "\t", row.names = F)
    
    ### store file path
    deg_filepath_vec <- c(file2write, deg_filepath_vec)
    chr_region_vec <- c(chr_region_tmp, chr_region_vec)
    aliquot_vec <- c(aliquot_tmp, aliquot_vec)
  }
}

