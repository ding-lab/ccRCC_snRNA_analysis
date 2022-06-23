# Yige Wu @WashU May 2022
## https://www.notion.so/Analyze-Previously-published-sc-snRNA-seq-data-292306b08ad54920adc5a0b36a186e0b

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
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "Matrix"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
options(future.globals.maxSize = 1000 * 1024^2)

# input ------------------------------------------------------
## input the barcode annotation
barcode_anno_df1 <- fread(data.table = F, input = "./Resources/Knowledge/Published_Data/Zhang_scRNA_PNAS_2021/tumor_anno.csv")
barcode_anno_df2 <- fread(data.table = F, input = "./Resources/Knowledge/Published_Data/Zhang_scRNA_PNAS_2021/normal_anno.csv")
barcode_anno_df3 <- fread(data.table = F, input = "./Resources/Knowledge/Published_Data/Zhang_scRNA_PNAS_2021/chRCC_anno.csv")

# combine barcode annotation ----------------------------------------------
barcode_anno_df <- rbind(barcode_anno_df1, barcode_anno_df2, barcode_anno_df3)
barcode_anno_df <- barcode_anno_df %>%
  mutate(SI_ID = str_split_fixed(string = cell, pattern = "_", n = 3)[,2]) %>%
  mutate(barcode = str_split_fixed(string = cell, pattern = "_", n = 3)[,3])

# process to seurat object ------------------------------------------------
files_input <- list.files(path = "./Resources/Knowledge/Published_Data/Zhang_scRNA_PNAS_2021/GSE159115_RAW/")
srat_list <- list()
for (filename_tmp in files_input) {
  count_mat_tmp <- Read10X_h5(filename = paste0("./Resources/Knowledge/Published_Data/Zhang_scRNA_PNAS_2021/GSE159115_RAW/", filename_tmp))
  si_id_tmp <- str_split_fixed(string = filename_tmp, pattern = "_", n = 8)[,3]
  ## get barcodes
  barcode_keep_df <- barcode_anno_df %>%
    filter(SI_ID == si_id_tmp) %>%
    filter(!(anno %in% c("ua", "Unknown")))
  barcodes_keep <- barcode_keep_df$barcode; length(barcodes_keep)
  ## create seurat object
  srat_tmp = CreateSeuratObject(counts = count_mat_tmp[,barcodes_keep], project= si_id_tmp, min.cells = 0)
  rm(count_mat_tmp)
  print(head(srat_tmp@meta.data))
  srat_tmp@meta.data$barcode <- rownames(srat_tmp@meta.data)
  metadata_df <- data.frame(srat_tmp@meta.data)
  metadata_df$row_id <- paste(1:nrow(metadata_df))
  metadata_df <- merge(x = metadata_df, 
                       y = barcode_keep_df, by.x = c("orig.ident", "barcode"), by.y = c("SI_ID", "barcode"), sort = F)
  rownames(metadata_df) <- metadata_df$barcode
  print(head(metadata_df))
  srat_tmp@meta.data <- metadata_df
  srat_list[[si_id_tmp]] <- srat_tmp
}
length(srat_list)
rm(srat_tmp)

# Run the standard workflow for visualization and clustering ------------
srat_list <- lapply(X = srat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
## integrate without anchor
srat_merged_obj <- merge(x = srat_list[[1]], y = srat_list[2:length(srat_list)], project = "merged")
rm(srat_list)

srat_merged_obj <- SCTransform(srat_merged_obj, vars.to.regress = c("nCount_RNA","pct_MT"), return.only.var.genes = F)
cat("Finished SCTransform!\n")
srat <- RunPCA(srat_merged_obj, npcs = 30, verbose = T)
cat("Finished RUNPCA!\n")
srat <- RunUMAP(srat_merged_obj, reduction = "pca", dims = 1:30)
cat("Finished RUNUMAP!\n")
srat <- FindNeighbors(srat_merged_obj, reduction = "pca", dims = 1:30, force.recalc = T)
cat("Finished FindNeighbors!\n")
srat <- FindClusters(srat_merged_obj, resolution = 0.5)
cat("Finished FindClusters!\n")

## save as RDS file
path_final_file <- paste0(dir_out, "Zhang_scRNA_PNAS_2021.Merged.", run_id, ".RDS")
saveRDS(object = srat_merged_obj, file = path_final_file, compress = T)
cat("Finished saving the output!\n")

