# Yige Wu @WashU Feb 2020
## fetch UMAP coordinates and cluster info for the tumor cell reclustering

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
## chromosome regions from which the CNVs are annotated here
chr_regions2process <- unique(ccrcc_cna_genes_df$chr_region)
chr_regions2process <- as.vector(chr_regions2process)
## set the directory with the DEG tables
dir_deg <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/findallmarkers/findallmarkers_expected_cnv_vs_not_in_malignant_nephron_epithelium/20200228.v1/"
deg_filepaths <- list.files(path = dir_deg, recursive = T) 
# loop for each chr_region ------------------------------------------------
deg_by_chr_region_df <- NULL
for (chr_region_tmp in chr_regions2process) {
  ## gather all the file paths
  chr_region_deg_filepaths <- deg_filepaths[grepl(pattern = paste0(chr_region_tmp, ".ExpectedCNV_vs_Neutral"), x = deg_filepaths)]
  chr_region_deg_filepaths
  
  ## input the files, filter by p.val.adj and bind
  for (filepath_tmp in chr_region_deg_filepaths) {
    deg_tmp <- fread(input = paste0(dir_deg, filepath_tmp), data.table = F)
    deg_tmp <- deg_tmp %>%
      filter(p_val < 0.05) %>%
      mutate(aliquot = str_split_fixed(string = filepath_tmp, pattern = "\\/", n = 2)[,1]) %>%
      mutate(chr_region = chr_region_tmp)
    deg_by_chr_region_df <- rbind(deg_tmp, deg_by_chr_region_df)
  }
}
unique(deg_by_chr_region_df$aliquot)
# write table -------------------------------------------------------------
write.table(x = deg_by_chr_region_df, file = paste0(dir_out, "FindMarkers.Wilcox.ExpectedCNV_vs_Neutral.", run_id, ".tsv"), sep = "\t", quote = F, row.names = F)