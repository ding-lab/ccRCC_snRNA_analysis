# Yige Wu @ WashU 20189 Aug
## for integrating the raw freature bc matrixs of multiple libraries together


###########################################
######## ANALYSIS
###########################################

# library -----------------------------------------------------------------
source("./Box/Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# Set path ----------------------------------------------------------------
matrix_dir <- "./Ding_Lab/Projects_Current/Multiple_Myeloma/MMY_Proteomics/Resources/scRNA/22933/22933-3_scRNA-1/raw_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

# Read input --------------------------------------------------------------
input <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(input) = barcode.names$V1
rownames(input) = feature.names$V1

bc_sums <- colSums(input)
bg_id <- names(bc_sums[bc_sums >= 300])

input = input[,bg_id]
kidney <- CreateSeuratObject(counts = input,project="Kidney",min.cells = 0)