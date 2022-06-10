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
# dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "doParallel"
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
## input the barcode to new cluster
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/annotate_clusterid_on_30ccRCCtumorcellreclustered_byresolution/20220406.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.20220406.v1.tsv")
## input clinical information
clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_case_clinical_data/20201125.v1/snRNA_ccRCC_Clinicl_Table.20201125.v1.tsv")
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input the gene set auto-correlation results
sigCorr_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/run_vision/getSignatureAutocorrelation_30ccRCC_tumorcellreclustered/20220411.v1/ccRCC.30ccRCC.TumorCellsReclustered.Vision.SignatureAutocorrelation.20220411.v1.tsv")

# preprocess --------------------------------------------------------------
## get barcodes to compare
clinical_filtered_df <- clinical_df %>% 
  filter(Case != "C3L-00359") %>%
  select(Case, Tumor_Stage_Pathological) 
cases_highstage <- clinical_filtered_df$Case[clinical_filtered_df$Tumor_Stage_Pathological %in% c("Stage III", "Stage IV")]
cases_lowstage <- clinical_filtered_df$Case[clinical_filtered_df$Tumor_Stage_Pathological %in% c("Stage I", "Stage II")]
alights_highstage <- metadata_df$Aliquot.snRNA[metadata_df$Case %in% cases_highstage & metadata_df$snRNA_available & metadata_df$Sample_Type == "Tumor"]; alights_highstage
alights_lowstage <- metadata_df$Aliquot.snRNA[metadata_df$Case %in% cases_lowstage & metadata_df$snRNA_available & metadata_df$Sample_Type == "Tumor"]; alights_lowstage
barcodes_highstage <- barcode2cluster_df$barcode[barcode2cluster_df$orig.ident %in% alights_highstage]
barcodes_lowstage <- barcode2cluster_df$barcode[barcode2cluster_df$orig.ident %in% alights_lowstage]

# prepare data to plot -----------------------------------------------------------------
# genesets_test <- colnames(sigScores)
genesets_test <- sigCorr_df$gene_set[sigCorr_df$FDR < 0.05]
registerDoParallel(cores = 20)
print("Start foreach!\n")
start_time <- Sys.time()
result_list<-foreach(g=genesets_test) %dopar% {
  sigscores_cluster1 <- sigScores[barcodes_highstage, g]
  sigscores_cluster2 <- sigScores[barcodes_lowstage, g]
  median_cluster1 <- median(sigscores_cluster1)
  median_diff <- median(sigscores_cluster1) - median(sigscores_cluster2)
  log2FC <- log2(median(sigscores_cluster1)/median(sigscores_cluster2))
  stat <- wilcox.test(x = sigscores_cluster1, y = sigscores_cluster2)
  p_val <- stat$p.value
  result_tmp <- list(c(p_val, median_diff, log2FC, median_cluster1))
  return(result_tmp)
}
print("Finish foreach!\n")
end_time <- Sys.time()
end_time - start_time 
results_df <- data.frame(matrix(data = unlist(result_list), ncol = 4, byrow = T))
print("Finish rbind.data.frame result_list!\n")
print(head(results_df))
colnames(results_df) <- c("p_val", "median_diff", "log2FC", "median_sigScore")
results_df$fdr <- p.adjust(p = results_df$p_val, method = "fdr")
results_df$gene_set <- genesets_test
print("Finish adding id columns!\n")

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "wilcox.vision_score.high_vs_low_stage.", run_id,".tsv")
write.table(x = results_df, file = file2write, quote = F, sep = "\t", row.names = F)




