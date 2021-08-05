# Yige Wu @WashU Oct 2019
## for running inferCNV using individual seruat objects

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

# set input & output directory ----------------------------------------------------
dir_infercnv <- "./Resources/snRNA_Processed_Data/InferCNV/"
dir_infercnv_inputs <- paste0(dir_infercnv, "inputs/")
dir.create(dir_infercnv_inputs)
dir_infercnv_counts <- paste0(dir_infercnv_inputs, "raw_counts_matrix/doublets_removed/")
dir.create(dir_infercnv_counts)

# input dependencies ------------------------------------------------
## input paths to the seurat objects
srat_paths_df <- fread(input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20210729.v1/Seurat_Object_Paths.20210729.v1.tsv", data.table = F)
srat_paths_df$Path_seurat_object
## input doublet prediction
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)

# input gene position file ------------------------------------------------
tab_gene_order_symbol = read.delim(file = "./Resources/snRNA_Processed_Data/InferCNV/inputs/gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt", header = FALSE,stringsAsFactors = FALSE)

# write raw counts matrix -------------------------------------------------
for (snRNA_aliquot_id_tmp in c("CPT0012280004", "CPT0012550012")) {
# for (snRNA_aliquot_id_tmp in unique(srat_paths_df$Aliquot.snRNA)) {
  path_infercnv_counts_out <- paste0(dir_infercnv_counts, snRNA_aliquot_id_tmp, ".RNA_Count.tsv")
  
  if (!file.exists(path_infercnv_counts_out)) {
    ## input the seurat object
    seurat_obj_path <- srat_paths_df$Path_katmai[srat_paths_df$Aliquot.snRNA == snRNA_aliquot_id_tmp]
    seurat_object <- readRDS(file = seurat_obj_path)
    
    ## get raw read count matrix
    raw_exp_mat <- seurat_object@assays$RNA@counts
    dim(raw_exp_mat)
    
    ## remove gene not mapped in the gene position file
    missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% tab_gene_order_symbol$V1)]
    missing_sym
    
    ## remove doublets
    barcodes_process <- colnames(raw_exp_mat)
    if (snRNA_aliquot_id_tmp %in% barcode2scrublet_df$Aliquot) {
      barcodes_doublets <- barcode2scrublet_df$Barcode[barcode2scrublet_df$Aliquot == snRNA_aliquot_id_tmp & barcode2scrublet_df$predicted_doublet]
      barcodes_process <- barcodes_process[!(barcodes_process %in% barcodes_doublets)]
    }
    
    ## only missing ~750 genes, romove them
    ## only keep the barcodes in the annotation files
    clean_exp_mat <- raw_exp_mat[!(rownames(raw_exp_mat) %in% missing_sym), barcodes_process]
    dim(clean_exp_mat)
    rm(raw_exp_mat)
    clean_exp_mat <- as.matrix(clean_exp_mat)
    write.table(x = clean_exp_mat, file = path_infercnv_counts_out, quote = F, row.names = T, col.names = T, sep = "\t")
  }
}

