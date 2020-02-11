#!/bin/bash

docker_image_name=sridnona/seurat_docker
dir_run=/diskmnt/Projects/ccRCC_scratch/
path_rscript=/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA_analysis/integration/integrate_30_seurat_objects/integrate_seurat_objects_on_cluster.R
dir_out=${dir_run}Resources/snRNA_Processed_Data/Analysis_Results/integration/30_aliquot_integration/run_integration/20200211.v1/

mkdir -p ${dir_out}
path_log_file=${dir_out}$0.$(date +%Y%m%d%H%M%S).log

docker run -v ${dir_run}:${dir_run} ${docker_image_name} Rscript ${path_rscript} &> ${path_log_file}&

echo "Log file is at "${path_log_file}
