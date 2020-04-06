# Yige Wu @WashU Oct 2019
## for running inferCNV using individual seruat objects

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set input & output directory ----------------------------------------------------
dir_infercnv <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/"
dir_infercnv_inputs <- paste0(dir_infercnv, "inputs/")
dir.create(dir_infercnv_inputs)
dir_infercnv_counts <- paste0(dir_infercnv_inputs, "raw_counts_matrix/")
dir.create(dir_infercnv_counts)
dir_infercnv_annotation <- paste0(dir_infercnv_inputs, "annotations_file/")
dir.create(dir_infercnv_annotation)

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200207.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(FACS == "") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object


# input gene position file ------------------------------------------------
tab_gene_order_symbol = read.delim(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/inputs/gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt", header = FALSE,stringsAsFactors = FALSE)

# write raw counts matrix -------------------------------------------------
for (snRNA_aliquot_id_tmp in unique(seurat_summary2process$Aliquot)) {
  path_infercnv_counts_out <- paste0(dir_infercnv_counts, snRNA_aliquot_id_tmp, ".RNA_Count.tsv")
  
  if (!file.exists(path_infercnv_counts_out)) {
    ## input the seurat object
    seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
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

