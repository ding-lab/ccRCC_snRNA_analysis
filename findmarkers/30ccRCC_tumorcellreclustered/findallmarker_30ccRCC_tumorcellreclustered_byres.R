# Yige Wu @WashU Mar 2022
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
  "data.table",
  "limma",
  "Seurat",
  "future",
  "future.apply",
  "doParallel",
  "ggplot2"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
# set up future for parallelization
plan("multiprocess", workers = 25)
options(future.globals.maxSize = 10000 * 1024^2)
# registerDoParallel(cores = 6)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out_parent <- makeOutDir_katmai(path_this_script)
dir_out <- paste0(dir_out_parent, run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_30_ccRCC_tumorcells/20220404.v1/ccRCC.34samples.Tumorcells.SeuratIntegrated.20220404.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")
logfc.threshold.run <- 0.25
min.pct.run <- 0.1
min.diff.pct.run <- 0
## input the barcode-to-cluster results
# barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/FindClusters_30_ccRCC_tumorcells_changeresolutions/20220405.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.20220405.v1.tsv")
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/FindClusters_30_ccRCC_tumorcells_changeresolutions/20220408.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.20220408.v1.tsv")

# process -----------------------------------------------------------------
resolutions_process <- colnames(barcode2cluster_df)[grepl(pattern = "integrated_snn_res", x = colnames(barcode2cluster_df))]
resolutions_process <- gsub(x = resolutions_process, pattern = "integrated_snn_res\\.", replacement = "")
result_df <- NULL
for (resolution_tmp in c("0.1", "0.5", "1", "2", "3", "4")) {
  path_markers <- paste0(dir_out_parent, "res.", resolution_tmp, ".tumorcellsreclustered.markers.logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
  if (file.exists(path_markers)) {
    markers <- fread(data.table = F, input = path_markers)
    cat(paste0("Markers for resolution ", resolution_tmp, "exists, reading!\n"))
  } else {
    cat(paste0("Markers for resolution ", resolution_tmp, "doesn't exist, running FindMarkers!\n"))
    
    srat@meta.data$cluster_test <- mapvalues(x = rownames(srat@meta.data), from = barcode2cluster_df$barcode, to = as.vector(barcode2cluster_df[, paste0("integrated_snn_res.", resolution_tmp)]))
    Idents(srat) <- "cluster_test"
    
    markers <- FindAllMarkers(object = srat, test.use = "wilcox", only.pos = T,
                              min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
    markers$gene_symbol <- rownames(markers)
    markers$resolution <- resolution_tmp
    write.table(x = markers, file = path_markers, quote = F, sep = "\t", row.names = F)
  }
  result_df <- rbind(result_df, markers)
}

# result_list<-foreach(resolution_tmp=c("0.1", "0.5", "1", "2", "3", "4")) %dopar% {
#   path_markers <- paste0(dir_out_parent, "res.", resolution_tmp, ".tumorcellsreclustered.markers.logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
#   if (file.exists(path_markers)) {
#     markers <- fread(data.table = F, input = path_markers)
#     cat(paste0("Markers for resolution ", resolution_tmp, "exists, reading!\n"))
#   } else {
#     cat(paste0("Markers for resolution ", resolution_tmp, "doesn't exist, running FindMarkers!\n"))
#     
#     srat@meta.data$cluster_test <- mapvalues(x = rownames(srat@meta.data), from = barcode2cluster_df$barcode, to = as.vector(barcode2cluster_df[, paste0("integrated_snn_res.", resolution_tmp)]))
#     Idents(srat) <- "cluster_test"
# 
#     markers <- FindAllMarkers(object = srat, test.use = "wilcox", only.pos = T,
#                            min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
#     markers$gene_symbol <- rownames(markers)
#     markers$resolution <- resolution_tmp
#     write.table(x = markers, file = path_markers, quote = F, sep = "\t", row.names = F)
#   }
#   return(markers)
# }
# result_df <- do.call(rbind.data.frame, result_list)

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "tumorcellsreclustered.markers.logfcthreshold.", 
                     logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run,
                     ".byresolution", ".tsv")
write.table(x = result_df, file = file2write, quote = F, sep = "\t", row.names = F)


