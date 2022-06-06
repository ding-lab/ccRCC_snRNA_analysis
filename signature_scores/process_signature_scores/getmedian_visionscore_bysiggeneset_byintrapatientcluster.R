# Yige Wu @WashU Mar 2021
## https://satijalab.org/seurat/archive/v3.0/integration.html
## also used references

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
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input ------------------------------------------------------
## input signature scores by barcode
sigScores <- readRDS("./Resources/Analysis_Results/signature_scores/run_vision/run_vision_on_30ccRCC_tumorcellreclustered/20220406.v1/ccRCC.30ccRCC.TumorCellsReclustered.Vision.sigScores.20220406.v1.RDS")
print("Finish reading the sigScores matrix!\n")
## input the gene set auto-correlation results
sigCorr_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/run_vision/getSignatureAutocorrelation_30ccRCC_tumorcellreclustered/20220411.v1/ccRCC.30ccRCC.TumorCellsReclustered.Vision.SignatureAutocorrelation.20220411.v1.tsv")
## input the manual intrapatient cluster
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)
## input the tumor cell intergrated barcode info
barcode_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/FindClusters_30_ccRCC_tumorcells_changeresolutions/20220405.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.20220405.v1.tsv")

# pre-process -----------------------------------------------------------------
genesets_test <- sigCorr_df$gene_set[sigCorr_df$FDR < 0.05]
barcode2tumorsubcluster_df <- barcode2tumorsubcluster_df %>%
  mutate(barcode_id_uniq = paste0(orig.ident, "_", barcode))
barcode_metadata_df <- barcode_metadata_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]) %>%
  mutate(barcode_id_uniq = paste0(orig.ident, "_", barcode_individual))
barcode_metadata_df <- merge(x = barcode_metadata_df, 
                             y = barcode2tumorsubcluster_df %>%
                               select(barcode_id_uniq, Cluster_Name) %>%
                               rename(intrapatient_cluster_name = Cluster_Name), 
                             by = c("barcode_id_uniq"), all.x = T)
which(is.na(barcode_metadata_df$intrapatient_cluster_name))

# process each tumor cluster in loop --------------------------------------
ic_name_tmp <- "C3L-00088-T1_C1"
id_names_process <- unique(barcode_metadata_df$intrapatient_cluster_name)
sigScores_med_df <- NULL
for (ic_name_tmp in id_names_process) {
  barcodes_tmp <- barcode_metadata_df$barcode[barcode_metadata_df$intrapatient_cluster_name == ic_name_tmp]
  sigScores_tmp <- sigScores[barcodes_tmp, genesets_test]
  sigScores_med_tmp <- apply(sigScores_tmp, 2, median)
  sigScores_med_df <- rbind(sigScores_med_df, t(sigScores_med_tmp))
}
sigScores_med_df <- as.data.frame(sigScores_med_df)
sigScores_med_df$intrapatient_cluster_name <- id_names_process

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "median_vision_score.sig_genesets.byintrapatient_tumor_cluster.", run_id,".tsv")
write.table(x = sigScores_med_df, file = file2write, quote = F, sep = "\t", row.names = F)


