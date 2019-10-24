#!/bin/bash

#docker_image_name=singlecellportal/infercnv
docker_image_name=sridnona/seurat_docker
dir_run=/diskmnt/Projects/ccRCC_scratch/
path_rscript=/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA_analysis/integration/integrate_seurat_objects_on_cluster.20191021.v1.R
dir_out=${dir_run}Resources/snRNA_Processed_Data/Analysis_Results/integration/integrate_seurat_objects_on_cluster/20191021.v1/

mkdir -p ${dir_out}
path_log_file=${dir_out}$0.$(date +%Y%m%d%H%M%S).log

docker run -v ${dir_run}:${dir_run} ${docker_image_name} Rscript ${path_rscript} &> ${path_log_file}&
