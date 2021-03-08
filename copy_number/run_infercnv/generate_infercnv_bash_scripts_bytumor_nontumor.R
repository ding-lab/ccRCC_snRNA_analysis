# Yige Wu @WashU Apr 2020
## for running inferCNV using individual seruat objects
## use immune and stromal cells as reference cells

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# set output directory------------------------------------------------------
dir_infercnv <- "./Resources/snRNA_Processed_Data/InferCNV/"
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
paths_srat <- fread(input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20200717.v1/Seurat_Object_Paths.20200717.v1.tsv", data.table = F)

## input barcode to cell type info
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200406.v1/30_aliquot_integration.barcode2celltype.20200406.v1.tsv", data.table = F)

snRNA_aliquot_id_tmp <- "CPT0023690004"
# write bash script to run docker image -----------------------------------
for (snRNA_aliquot_id_tmp in unique(paths_srat$Aliquot)) {
  ## make the filename for the bash script
  path_infercnv_script <- paste0(dir_infercnv_scripts, "run_sample.", snRNA_aliquot_id_tmp, ".sh")
  
  if (!file.exists(path_infercnv_script)) {
    ## input annotation file
    path_annotations_file_out <- paste0(dir_infercnv_annotation_out, snRNA_aliquot_id_tmp, ".Barcode_Annotation.txt")
    anno_tab <- fread(input = path_annotations_file_out, data.table = F, col.names = c("barcode", "infercnv_group"))
    ref_infercnv_group <- "Ref"
    ## start writing to the bash script
    sink(path_infercnv_script)
    cat("#!/bin/bash\n")
    cat("\n")
    cat(paste0("snRNA_aliquot_id=", snRNA_aliquot_id_tmp, "\n"))
    cat("analysis_mode=subclusters\n")
    cat("dir_run=/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/\n")
    cat("docker_image_name=singlecellportal/infercnv\n")
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
    cat(paste0("ref_group_names=",paste0(ref_infercnv_group, collapse = ","), "\n"))
    # cat(paste0("num_ref_groups=", length(ref_celltype_shorter), "\n"))
    cat(paste0("docker run -v ${dir_run}:${dir_run} --user $(id -u):$(id -g) --ulimit stack=8277716992:8277716992 ${docker_image_name} inferCNV.R \\\
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





