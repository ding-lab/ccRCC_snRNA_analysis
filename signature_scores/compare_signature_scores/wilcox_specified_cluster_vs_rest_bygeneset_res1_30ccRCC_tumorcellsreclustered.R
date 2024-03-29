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

# input ------------------------------------------------------
## input signature scores by barcode
sigScores <- readRDS("./Resources/Analysis_Results/signature_scores/run_vision/run_vision_on_30ccRCC_tumorcellreclustered/20220406.v1/ccRCC.30ccRCC.TumorCellsReclustered.Vision.sigScores.20220406.v1.RDS")
print("Finish reading the sigScores matrix!\n")
## input the barcode to new cluster
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/annotate_clusterid_on_30ccRCCtumorcellreclustered_byresolution/20220406.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.20220406.v1.tsv")

# prepare data to plot -----------------------------------------------------------------
barcode2cluster_df$clusterid_test <- barcode2cluster_df$integrated_snn_res.1
cluster1_process <- c(0, 3, 7, 12)
genesets_test <- colnames(sigScores)
registerDoParallel(cores = 20)

cat(paste0("Start compare, ", paste0(cluster1_process, collapse = ","), " vs. ", "the rest", "!\n"))

barcodes_cluster1 <- barcode2cluster_df$barcode[barcode2cluster_df$clusterid_test %in% cluster1_process]
barcodes_cluster2 <- barcode2cluster_df$barcode[!(barcode2cluster_df$clusterid_test %in% cluster1_process)]

print("Start foreach!\n")
start_time <- Sys.time()
result_list<-foreach(g=genesets_test) %dopar% {
  sigscores_cluster1 <- sigScores[barcodes_cluster1, g]
  sigscores_cluster2 <- sigScores[barcodes_cluster2, g]
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
result_tmp_df <- data.frame(matrix(data = unlist(result_list), ncol = 4, byrow = T))
print("Finish rbind.data.frame result_list!\n")
print(head(result_tmp_df))

colnames(result_tmp_df) <- c("p_val", "median_diff", "log2FC", "median_sigScore")
result_tmp_df$fdr <- p.adjust(p = result_tmp_df$p_val, method = "fdr")
result_tmp_df$gene_set <- genesets_test
print("Finish adding id columns!\n")

# save output -------------------------------------------------------------
## set run id
run_id <- "0_3_7_12_vs_others.v1"
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "wilcox.res1.30ccRCCtumorcellreclustered,", run_id,".tsv")
write.table(x = result_tmp_df, file = file2write, quote = F, sep = "\t", row.names = F)




