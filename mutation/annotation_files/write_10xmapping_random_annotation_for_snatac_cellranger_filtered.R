# Yige Wu @WashU Dec 2020
## for generating the barcode annotation files (assigned to one random group) for 10Xmapping
## only use barcodes of after QC barcodes so that 10Xmapping won't run forever

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
## input the barcode info
umap_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Barcode_Annotation/UMAP/UMAP_data_13_snATAC_Normal_epithelialCells_reclustered.20201209.tsv")
## input meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# input seurat object -----------------------------------------------------
for (aliquot in unique(umap_df$dataset)) {
  easy_id <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquot]
  anno_tab_tmp <- umap_df %>%
    filter(dataset == aliquot) %>%
    mutate(barcode = str_split_fixed(string = V1, pattern = "_", n = 2)[,2]) %>%
    select(barcode) %>%
    mutate(random_group = "0")
  write.table(x = anno_tab_tmp, file = paste0(dir_out, aliquot, "_snATAC", "_AfterQC_Barcodes.tsv"), quote = F, row.names = F, sep = "\t", col.names = F)
}


