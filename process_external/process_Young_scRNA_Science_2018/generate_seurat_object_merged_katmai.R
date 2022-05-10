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
matrix.path <- "./Resources/Knowledge/Published_Data/Young_scRNA_Science_2018/aat1699_DataS1/tableOfCounts.mtx"
features.path <- "./Resources/Knowledge/Published_Data/Young_scRNA_Science_2018/aat1699_DataS1/tableOfCounts_rowLabels.tsv"
barcode.path <- "./Resources/Knowledge/Published_Data/Young_scRNA_Science_2018/aat1699_DataS1/tableOfCounts_colLabels.tsv"
input <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = T,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = T,stringsAsFactors = FALSE)
## input the barcode annotation
barcode_anno_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Published_Data/Young_scRNA_Science_2018/aat1699-young-tabless1-s12-revision2.xlsx", 
                                     sheet = "TableS11 - Cell manifest", skip = 1)

# process to seurat object ------------------------------------------------
## check dimensions and duplicates in dimension names
dim(input); nrow(barcode.names)
which(duplicated(barcode.names$Barcode)) ## the "Barcode" column has duplicates
which(duplicated(barcode.names$DropletID)) ## the "DropletID" column doesn't have duplicates
which(duplicated(feature.names$Symbol)) ## the "Symbol" column has duplicates
feature.names[which(duplicated(feature.names$Symbol)),]
## add dimension names
colnames(input) = barcode.names$DropletID
rownames(input) = feature.names$GeneLabel
## only get QC-passed droplet ids
dropletids_pass <- barcode_anno_df$DropletID[barcode_anno_df$QCpass]
length(dropletids_pass) # 72501
## create seurat object
dim(input[,dropletids_pass])
srat = CreateSeuratObject(counts = input[,dropletids_pass], project= "Young_scRNA_Science_2018",min.cells = 0)
rm(input)
## add barcode meta data
barcode_anno_df <- data.frame(barcode_anno_df)
srat@meta.data$DropletID <- rownames(srat@meta.data)
metadata_df <- data.frame(srat@meta.data)
metadata_df$row_id <- 1:nrow(metadata_df)
metadata_df <- merge(x = metadata_df, 
                     y = barcode_anno_df %>%
                       filter(QCpass), by = c("DropletID"), sort = F)
nrow(metadata_df)
rownames(metadata_df) <- metadata_df$DropletID
srat@meta.data <- metadata_df
View(srat@meta.data)
## add feature meta.data
feature.names <- feature.names %>%
  mutate(feature_name = gsub(pattern = "_ENSG", replacement = "-ENSG", x = GeneLabel))
rownames(feature.names) <- feature.names$feature_name
srat@assays$RNA@meta.features <- feature.names[rownames(srat@assays$RNA@counts),]


# rest of the pipeline ----------------------------------------------------

## scale data with all the features
path_final_file <- paste0(dir_out, "Young_scRNA_Science_2018.Merged.", run_id, ".RDS")
path_sct_file <- paste0(dir_out, "Young_scRNA_Science_2018.SCT." , ".RDS")

if (!file.exists(path_sct_file)) {
  ## rest of the pipeline
  srat <- SCTransform(srat, vars.to.regress = c("nCount_RNA","MTfrac"), return.only.var.genes = F)
  saveRDS(object = srat, file = path_sct_file, compress = T)
  cat("Finished Writing SCTransform!\n")
  
  srat <- RunPCA(srat, npcs = 30, verbose = FALSE)
  cat("Finished RUNPCA!\n")
  srat <- RunUMAP(srat, reduction = "pca", dims = 1:30)
  cat("Finished RUNUMAP!\n")
  srat <- FindNeighbors(srat, reduction = "pca", dims = 1:30, force.recalc = T)
  cat("Finished FindNeighbors!\n")
  srat <- FindClusters(srat, resolution = 0.5)
  cat("Finished FindClusters!\n")
  
} else {
  srat <- readRDS(file = path_sct_file)
}
## save as RDS file
file2write <- paste0(dir_out, "Young_scRNA_Science_2018.Seurat.RDS")
saveRDS(object = srat, file = file2write, compress = T)
cat("Finished saving the output!\n")


