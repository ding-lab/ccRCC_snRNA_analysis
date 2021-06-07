# Yige Wu @ WashU 2021 Jun
## annotate sample copy number profile (3p, 5q, 14q)
## CNV: https://wustl.box.com/s/vlde6w791k81q1aibvgs0a487zxleij6
# Purity and ploidy: https://wustl.box.com/s/jel5krgvnvlq5z32vdg4itdcq6znuqzn
# From UMich

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
# library(future)
# plan("multiprocess", workers = 4)
library(foreach)
library(doParallel)
doParallel::registerDoParallel(cl = makeCluster(4), cores = 5)
getDoParWorkers()
options(future.globals.maxSize = 5 * 1024^3) # for 5 Gb RAM

# input dependencies ------------------------------------------------------
## input CNA matrix
cna_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/WGS_CNV_Somatic/Absolute_cnv/c3-ccrcc-combined-cnvex-lr_v1.0.csv")
## input snRNA sample set
metadata_df <- fread("./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv", data.table = F)
## input peaks
peaks_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Peak_Annotation/All_peaks_annotated_26snATAC_merged_obj.20210607.tsv")
## input the barcode info
barcodes_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Barcode_Annotation/Cell_type_annotation_snATAC.20210604.tsv")

# preprocess ----------------------------
## get barcodes to process
barcodes_df <- barcodes_df %>%
  mutate(id_bc = paste0(Sample, "_", Barcode))
barcodes_df$Case <- mapvalues(x = barcodes_df$Sample, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Case))
barcodes_df$Sample_type <- mapvalues(x = barcodes_df$Sample, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Sample_Type))
# table(barcodes_df$Cell_type)
# barcodes_df %>%
#   filter(Sample_type == "Tumor" & Cell_type == "PT")
barcodes_process <- barcodes_df$id_bc[barcodes_df$Cell_type %in% c("Tumor", "PT")]
cases_process <- barcodes_df$Case[barcodes_df$Cell_type %in% c("Tumor", "PT")]
rm(barcodes_df)
## get genes to process
peaks_df <- peaks_df %>%
  rename(Gene = SYMBOL)
genes_process <- unique(peaks_df$Gene)
genes_process <- cna_df$gene_name[cna_df$gene_name %in% genes_process]

# preprocess mean CNV values per gene--------------------------------------------------------------
start_time <- Sys.time()
foreach (gene_tmp = genes_process[1:50]) %dopar% {
  ## filter the CNVs
  cna_filtered_df <- cna_df[cna_df$gene_name == gene_tmp,]
  cna_filtered_df <- cna_filtered_df[!duplicated(cna_filtered_df$gene_name),]
  cna_filtered_vec <- unlist(cna_filtered_df[, unique(cases_process)])
  cna_bybc_vec <- cna_filtered_vec[cases_process]
  cna_bybc_vec[celltypes_process == "PT"] <- 0
  cna_t_df <- data.frame(cnv_gene = cna_bybc_vec)
  rownames(cna_t_df) <- barcodes_process
  # cna_t_df[1:5,]
  
  # write outupt -------------------------------------------------------------
  ## based on this, I'll save as RDS compressed: https://waterdata.usgs.gov/blog/formats/
  file2write <- paste0(dir_out, gene_tmp, ".CNV.ByATACBarcode",  ".RDS")
  saveRDS(object = cna_t_df, file = file2write, compress = T)
  # file2write <- paste0(dir_out, gene_tmp, ".CNV.ByATACBarcode",  ".tsv")
  # write.table(x = cna_t_df, file = file2write, quote = F, sep = "\t", row.names = T)
}
end_time <- Sys.time()
end_time - start_time