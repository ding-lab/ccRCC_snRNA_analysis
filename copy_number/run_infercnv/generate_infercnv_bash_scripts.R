# Yige Wu @WashU Apr 2020
## for running inferCNV using individual seruat objects
## use immune and stromal cells as reference cells

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
run_id <- "20200408.v1"

# set output directory------------------------------------------------------
dir_infercnv <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/"
dir_infercnv_inputs <- paste0(dir_infercnv, "inputs/")
dir.create(dir_infercnv_inputs)
dir_infercnv_counts <- paste0(dir_infercnv_inputs, "raw_counts_matrix/")
dir.create(dir_infercnv_counts)

dir_infercnv_scripts <- paste0(dir_infercnv, "scripts", "/", "testing/", "docker_call/","Individual.", run_id, "/")
dir.create(dir_infercnv_scripts)

dir_infercnv_annotation <- paste0(dir_infercnv_inputs, "annotations_file/")
dir.create(dir_infercnv_annotation)
dir_infercnv_annotation_out <- paste0(dir_infercnv_annotation, "Individual.", run_id, "/")
dir.create(dir_infercnv_annotation_out)

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
## input barcode to cell type info
barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200406.v1/30_aliquot_integration.barcode2celltype.20200406.v1.tsv", data.table = F)

snRNA_aliquot_id_tmp <- "CPT0023690004"
# write bash script to run docker image -----------------------------------
for (snRNA_aliquot_id_tmp in unique(seurat_summary2process$Aliquot)) {
  ## make the filename for the bash script
  path_infercnv_script <- paste0(dir_infercnv_scripts, "run_sample.", snRNA_aliquot_id_tmp, ".sh")
  
  if (!file.exists(path_infercnv_script)) {
    ## input annotation file
    path_annotations_file_out <- paste0(dir_infercnv_annotation_out, snRNA_aliquot_id_tmp, ".Barcode_Annotation.txt")
    
    ref_barcode2celltype_tmp <- barcode2celltype_df %>%
      filter(orig.ident == snRNA_aliquot_id_tmp) %>%
      filter(!(Cell_type.shorter %in% c("Tumor cells", "Unknown")))
    ref_celltype_shorter <- unique(ref_barcode2celltype_tmp$Cell_type.shorter)
    ## start writing to the bash script
    sink(path_infercnv_script)
    cat("#!/bin/bash\n")
    cat("\n")
    cat(paste0("snRNA_aliquot_id=", snRNA_aliquot_id_tmp, "\n"))
    cat("analysis_mode=subclusters\n")
    cat("dir_run=/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/\n")
    # cat("docker_image_name=singlecellportal/infercnv\n")
    cat("docker_image_name=trinityctat/infercnv:latest\n")
    
    trinityctat/infercnv:latest
    cat(paste0("dir_infercnv=${dir_run}Resources/snRNA_Processed_Data/InferCNV/\n"))
    cat(paste0("dir_infercnv_outputs=${dir_infercnv}outputs/\n"))
    cat(paste0("dir_infercnv_inputs=${dir_infercnv}inputs/\n"))
    cat(paste0("run_id=", run_id, "\n"))
    cat(paste0("dir_infercnv_outputs_by_run=${dir_infercnv_outputs}Individual.${run_id}/\n"))
    cat("mkdir -p ${dir_infercnv_outputs_by_run}\n")
    cat(paste0("dir_output=${dir_infercnv_outputs_by_run}${snRNA_aliquot_id}/", "\n"))
    cat("mkdir -p ${dir_output}\n")
    
    cat(paste0("dir_raw_count_file=${dir_infercnv_inputs}raw_counts_matrix/", "\n"))
    cat(paste0("path_raw_count_file=${dir_raw_count_file}${snRNA_aliquot_id}.RNA_Count.tsv", "\n"))
    cat(paste0("dir_annotation_file=${dir_infercnv_inputs}annotations_file/Individual.${run_id}/", "\n"))
    cat(paste0("path_annotation_file=${dir_annotation_file}${snRNA_aliquot_id}.Barcode_Annotation.txt", "\n"))
    cat(paste0("path_log_file=${dir_output}${snRNA_aliquot_id}.$(date +%Y%m%d%H%M%S).log", "\n"))
    cat(paste0("path_gene_order_file=${dir_infercnv_inputs}gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt", "\n"))
    cat(paste0("cutoff=0.04", "\n"))
    cat(paste0("ref_group_names=","Ref", "\n"))
    cat(paste0("num_ref_groups=", length(ref_celltype_shorter), "\n"))
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
        --num_ref_groups=${num_ref_groups} \\\
        --ref_group_names=${ref_group_names} &> ${path_log_file}&", "\n"))
    cat(paste0("\n"))
    cat(paste0("\n"))
    cat(paste0("\n"))
    cat(paste0("\n"))
    cat(paste0("\n"))
    sink()
  }
}





