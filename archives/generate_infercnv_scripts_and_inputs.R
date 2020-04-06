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

# set run id ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set integration id and output directory------------------------------------------------------
integration_id <- "20191021.v1"
dir_infercnv_annotation_out <- paste0(dir_infercnv_annotation, "integration.", integration_id, "/")
dir.create(dir_infercnv_annotation_out)

# input cell type assignment table ----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191022.v1.tsv", data.table = F)

# get ref groups--------------------------------------------------------
tumor_group_names <- cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Malignant == "Yes"]
tumor_group_names <- as.character(tumor_group_names)
tumor_group_names

ref_group_names <- cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Malignant == "No"]
ref_group_names
## [1]  2  5  8  9 10 11 13 14 16 17 18 19 20 21

# input integrated seurat object ------------------------------------------
int_seurat_object <- readRDS("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal_Integrated.20191021.v1.RDS")
int_seurat_object@meta.data$orig.ident %>% unique()
# [1] "CPT0075140002" "CPT0001260013" "CPT0086350004" "CPT0001180011" "CPT0025890002" "CPT0010110013" "CPT0019130004" "CPT0020120013" "CPT0001220012"
# [10] "CPT0014450005"

# write raw counts matrix -------------------------------------------------
for (snRNA_aliquot_id in unique(int_seurat_object@meta.data$orig.ident)) {
  path_infercnv_counts_out <- paste0(dir_infercnv_counts, snRNA_aliquot_id, ".RNA_Count.tsv")
  
  if (!file.exists(path_infercnv_counts_out)) {
    # subset integrated object by sample ------------------------------------------------------
    seurat_object <- Seurat::SubsetData(int_seurat_object, idents = snRNA_aliquot_id)
    
    # get cells 2 process -----------------------------------------------------
    sample_barcodes_int <- rownames(seurat_object@meta.data[seurat_object@meta.data$orig.ident == snRNA_aliquot_id,])
    sample_barcodes_int
    
    # get raw read count matrix -----------------------------------------------
    raw_exp_mat <- seurat_object@assays$RNA@counts[, sample_barcodes_int]
    dim(raw_exp_mat)
    
    # change the barcode to exclude the underline -----------------------------
    sample_barcodes_sim <- as.vector(str_split_fixed(string = sample_barcodes_int, pattern = "_", n = 2)[,1])
    colnames(raw_exp_mat) <- sample_barcodes_sim
    
    # expression data - remove gene not mapped in the gene position file
    tab_gene_order_symbol = read.delim(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/inputs/gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt", header = FALSE,stringsAsFactors = FALSE)
    missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% tab_gene_order_symbol$V1)]
    missing_sym
    
    ## only missing ~750 genes, romove them
    ## only keep the barcodes in the annotation files
    clean_exp_mat <- raw_exp_mat[!(rownames(raw_exp_mat) %in% missing_sym), ]
    rm(raw_exp_mat)
    clean_exp_mat <- as.matrix(clean_exp_mat)
    write.table(x = clean_exp_mat, file = path_infercnv_counts_out, quote = F, row.names = T, col.names = T, sep = "\t")
    
    # check if the matrix is good for running ---------------------------------
    rownames(clean_exp_mat) %in% tab_gene_order_symbol$V1 %>% all # Check if only gene with position left. Should be TRUE
    rownames(clean_exp_mat) %>% duplicated %>% any # Check if duplicate still exist. Should be FALSE
  }
}

# write annotation files -------------------------------------------------
for (snRNA_aliquot_id in unique(int_seurat_object@meta.data$orig.ident)) {
  path_annotations_file_out <- paste0(dir_infercnv_annotation_out, snRNA_aliquot_id, ".Barcode_Annotation.txt")
  
  if (!file.exists(path_annotations_file_out)) {
    seurat_object <- Seurat::SubsetData(int_seurat_object, idents = snRNA_aliquot_id)
    
    anno_tab <- seurat_object@meta.data[seurat_object@meta.data$orig.ident == snRNA_aliquot_id,]
    anno_tab$barcode_int <- rownames(anno_tab)
    anno_tab <- anno_tab %>%
      mutate(barcode = str_split_fixed(string = barcode_int, pattern = "_", n = 2)[,1]) %>%
      select(barcode, seurat_clusters)
    nrow(anno_tab)
    write.table(x = anno_tab, file = path_annotations_file_out, quote = F, sep = "\t", row.names = F, col.names = F)
  }
}

dir_infercnv_scripts <- paste0(dir_infercnv, "InferCNV", "/", "testing/", "docker_call/","integration.", integration_id, "/")
dir.create(dir_infercnv_scripts)
# write bash script to run docker image -----------------------------------
for (snRNA_aliquot_id in unique(int_seurat_object@meta.data$orig.ident)) {
  path_infercnv_script <- paste0(dir_infercnv_scripts, "run_sample.", snRNA_aliquot_id, ".sh")
  
  if (!file.exists(path_infercnv_script)) {
    # subset integrated object by sample ------------------------------------------------------
    seurat_object <- Seurat::SubsetData(int_seurat_object, idents = snRNA_aliquot_id)
    
    ref_group_names2p <- intersect(ref_group_names, as.character(unique(seurat_object@meta.data$seurat_clusters)))
    
    sink(path_infercnv_script)
    cat("#!/bin/bash\n")
    cat("\n")
    cat(paste0("snRNA_aliquot_id=", snRNA_aliquot_id, "\n"))
    cat("analysis_mode=subclusters\n")
    cat("dir_run=/diskmnt/Projects/ccRCC_scratch/\n")
    cat("docker_image_name=singlecellportal/infercnv\n")
    cat(paste0("dir_infercnv=${dir_run}Resources/snRNA_Processed_Data/InferCNV/\n"))
    cat(paste0("dir_infercnv_outputs=${dir_infercnv}outputs/\n"))
    cat(paste0("dir_infercnv_inputs=${dir_infercnv}inputs/\n"))
    cat(paste0("integration_id=", integration_id, "\n"))
    cat(paste0("dir_infercnv_outputs_by_int=${dir_infercnv_outputs}integration.${integration_id}/\n"))
    cat("mkdir -p ${dir_infercnv_outputs_by_int}\n")
    cat(paste0("dir_output=${dir_infercnv_outputs_by_int}${snRNA_aliquot_id}/", "\n"))
    cat("mkdir -p ${dir_output}\n")
    
    cat(paste0("dir_raw_count_file=${dir_infercnv_inputs}raw_counts_matrix/", "\n"))
    cat(paste0("path_raw_count_file=${dir_raw_count_file}${snRNA_aliquot_id}.RNA_Count.tsv", "\n"))
    cat(paste0("dir_annotation_file=${dir_infercnv_inputs}annotations_file/integration.${integration_id}/", "\n"))
    cat(paste0("path_annotation_file=${dir_annotation_file}${snRNA_aliquot_id}.Barcode_Annotation.txt", "\n"))
    cat(paste0("path_log_file=${dir_output}${snRNA_aliquot_id}.$(date +%Y%m%d%H%M%S).log", "\n"))
    cat(paste0("path_gene_order_file=/diskmnt/Projects/ccRCC_scratch/Resources/snRNA_Processed_Data/InferCNV/inputs/gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt", "\n"))
    cat(paste0("cutoff=0.04", "\n"))
    cat(paste0("ref_group_names=", paste0(ref_group_names, collapse = ","), "\n"))
    cat(paste0("docker run -v ${dir_run}:${dir_run} ${docker_image_name} inferCNV.R \\\
        --analysis_mode=${analysis_mode} \\\
        --raw_counts_matrix=${path_raw_count_file} \\\
        --annotations_file=${path_annotation_file} \\\
        --gene_order_file=${path_gene_order_file} \\\
        --cutoff=${cutoff} \\\
        --out_dir=${dir_output} \\\
        --cluster_by_groups \\\
        --denoise \\\
        --HMM \\\
        --num_threads=20 \\\
        --ref_group_names=${ref_group_names} &> ${path_log_file}&", "\n"))
    cat(paste0("\n"))
    cat(paste0("\n"))
    cat(paste0("\n"))
    cat(paste0("\n"))
    cat(paste0("\n"))
    sink()
  }
}





