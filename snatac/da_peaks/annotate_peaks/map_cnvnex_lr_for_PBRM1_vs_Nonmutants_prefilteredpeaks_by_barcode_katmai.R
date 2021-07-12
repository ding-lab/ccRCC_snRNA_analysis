# Yige Wu @ WashU 2021 Jul
## annotate sample copy number profile for each peak
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
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 5 * 1024^3) # for 5 Gb RAM

# input dependencies ------------------------------------------------------
## input CNA matrix
cna_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/WGS_CNV_Somatic/Absolute_cnv/c3-ccrcc-combined-cnvex-lr_v1.0.csv")
## input snRNA sample set
metadata_df <- fread("./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv", data.table = F)
## input peaks
peaks_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/PBRM1_Specific/PBRM1_vsNon_mutants.Filtered_peaks_byMinPct.0.1.20210712.tsv")
peak2gene_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Peak_Annotation/28_snATACmerged_allPeaks.Annotated.20210712.tsv")
## input the barcode info
barcodes_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Barcode_Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv")

# preprocess ----------------------------
## get barcodes to process
barcodes_df <- barcodes_df %>%
  rename(id_bc = V1) %>%
  rename(Cell_type = cell_type) %>%
  rename(Sample = dataset)
  # mutate(id_bc = paste0(Sample, "_", Barcode))
barcodes_df$Case <- mapvalues(x = barcodes_df$Sample, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Case))
table(barcodes_df$Case)
barcodes_df$Sample_type <- mapvalues(x = barcodes_df$Sample, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Sample_Type))
table(barcodes_df$Sample_type)
# barcodes_df %>%
#   filter(Sample_type == "Tumor" & Cell_type == "PT")
idx_process <- (barcodes_df$Sample_type == "Tumor" & barcodes_df$Cell_type == "Tumor")
barcodes_process <- barcodes_df$id_bc[idx_process]
cases_process <- barcodes_df$Case[idx_process]
celltypes_process <- barcodes_df$Cell_type[idx_process]
## get genes to process
peak2gene_filtered_df <- peak2gene_df %>%
  filter(peak %in% peaks_df$peak) %>%
  rename(SYMBOL = Gene)
genes_process <- unique(peak2gene_filtered_df$SYMBOL)
genes_process <- genes_process[genes_process %in% cna_df$gene_name]

# preprocess mean CNV values per gene--------------------------------------------------------------
## filter the CNVs
cna_filtered_df <- cna_df[cna_df$gene_name %in% genes_process,]
## preprocess the CNV data frame
colnames_old <- colnames(cna_filtered_df)
colnames_new <- str_split_fixed(string = colnames_old, pattern = "\\.", n = 4)[,1]
colnames(cna_filtered_df) <- colnames_new
cna_filtered_df <- cna_filtered_df[!duplicated(cna_filtered_df$gene_name),]
## add peak
cna_bypeak_df <- merge(x = cna_filtered_df, 
                       y = peak2gene_filtered_df %>%
                         select(peak, SYMBOL),
                       by.x = c("gene_name"), by.y = c("SYMBOL"), all.x = T)

cna_bybc_bypeak_df <- cna_bypeak_df[, cases_process]
rownames(cna_bybc_bypeak_df) <- cna_bypeak_df$peak
colnames(cna_bybc_bypeak_df) <- barcodes_process
cna_bybc_bypeak_df[1:5, 1:5]
cna_t_mat <- t(as.matrix(cna_bybc_bypeak_df))
cna_t_mat[1:5, 1:5]

# write table -------------------------------------------------------------
file2write <- paste0(dir_out, "Barcode2PBRM1_vs_Nonmutants_PrefilteredPeak.CNV.", run_id, ".RDS")
saveRDS(object = cna_t_mat, file = file2write, compress = T)
