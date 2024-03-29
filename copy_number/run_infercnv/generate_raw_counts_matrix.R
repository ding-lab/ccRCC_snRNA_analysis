# Yige Wu @WashU Oct 2019
## for running inferCNV using individual seruat objects

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(dplyr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# set input & output directory ----------------------------------------------------
dir_infercnv <- "./Resources/snRNA_Processed_Data/InferCNV/"
dir_infercnv_inputs <- paste0(dir_infercnv, "inputs/")
dir.create(dir_infercnv_inputs)
dir_infercnv_counts <- paste0(dir_infercnv_inputs, "raw_counts_matrix/")
dir.create(dir_infercnv_counts)
dir_infercnv_annotation <- paste0(dir_infercnv_inputs, "annotations_file/")
dir.create(dir_infercnv_annotation)

# input seurat processing summary ------------------------------------------------
srat_paths_df <- fread(input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20210729.v1/Seurat_Object_Paths.20210729.v1.tsv", data.table = F)
srat_paths_df$Path_seurat_object

# input gene position file ------------------------------------------------
tab_gene_order_symbol = read.delim(file = "./Resources/snRNA_Processed_Data/InferCNV/inputs/gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt", header = FALSE,stringsAsFactors = FALSE)

# write raw counts matrix -------------------------------------------------
for (snRNA_aliquot_id_tmp in unique(srat_paths_df$Aliquot)) {
  path_infercnv_counts_out <- paste0(dir_infercnv_counts, snRNA_aliquot_id_tmp, ".RNA_Count.tsv")
  
  if (!file.exists(path_infercnv_counts_out)) {
    ## input the seurat object
    seurat_obj_path <- srat_paths_df$Path_seurat_object[srat_paths_df$Aliquot == snRNA_aliquot_id_tmp]
    seurat_object <- readRDS(file = seurat_obj_path)
    
    ## get raw read count matrix
    raw_exp_mat <- seurat_object@assays$RNA@counts
    dim(raw_exp_mat)
    
    ## remove gene not mapped in the gene position file
    missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% tab_gene_order_symbol$V1)]
    missing_sym
    
    ## only missing ~750 genes, romove them
    ## only keep the barcodes in the annotation files
    clean_exp_mat <- raw_exp_mat[!(rownames(raw_exp_mat) %in% missing_sym), ]
    rm(raw_exp_mat)
    clean_exp_mat <- as.matrix(clean_exp_mat)
    write.table(x = clean_exp_mat, file = path_infercnv_counts_out, quote = F, row.names = T, col.names = T, sep = "\t")
  }
}

