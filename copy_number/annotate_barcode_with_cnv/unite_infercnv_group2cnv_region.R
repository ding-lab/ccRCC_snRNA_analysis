# Yige Wu @WashU May 2020

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
## set infercnv output directory
dir_infercnv_output <- "./Resources/snRNA_Processed_Data/InferCNV/outputs/"
infercnv_run_id <- "Individual.20200305.v1"
dir_infercnv_run <- paste0(dir_infercnv_output, infercnv_run_id, "/")
## input id meta data 
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")

# read file by aliquot ----------------------------------------------------
## get aliquots to process
aliquots2process <- list.files(dir_infercnv_run)
aliquots2process <- aliquots2process[grepl(pattern = "CPT", x = aliquots2process)]
group2cnv_region_df <- NULL
for (aliquot_tmp in aliquots2process) {
  ## input infercnv cell grouping
  file2read <- paste0(dir_infercnv_run, aliquot_tmp, "/", "HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat")
  barcode2group_tmp <- fread(data.table = F, input = file2read)
  barcode2group_tmp$Id_Aliquot <- aliquot_tmp
  group2cnv_region_df <- rbind(group2cnv_region_df, barcode2group_tmp)
}
nrow(group2cnv_region_df)
group2cnv_region_df$Id_Aliquot_WU <- mapvalues(x = group2cnv_region_df$Id_Aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "InferCNV_Group2CNV_Region.", run_id, ".tsv")
write.table(x = group2cnv_region_df, file = file2write, sep = "\t", row.names = F, quote = F)

