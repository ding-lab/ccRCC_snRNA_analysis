#!/bin/bash

#docker_image_name=singlecellportal/infercnv
docker_image_name=sridnona/seurat_docker
dir_run=/diskmnt/Projects/ccRCC_scratch/
dir_rscript=/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA_analysis/recluster/
name_rscript=process_immune_cells_on_cluster.20191022.v1.R
path_rscript=${dir_rscript}${name_rscript}
dir_resources=${dir_run}Resources/
dir_analysis_results=${dir_resources}Analysis_Results/
dir_recluster_results=${dir_analysis_results}recluster/

dir_out_parent=${dir_recluster_results}process_immune_cells_on_cluster/
mkdir -p ${dir_out_parent}
dir_out=${dir_out_parent}20191022.v1/
mkdir -p ${dir_out}

path_log_file=${dir_out}$0.$(date +%Y%m%d%H%M%S).log

docker run -v ${dir_run}:${dir_run} ${docker_image_name} Rscript ${path_rscript} &> ${path_log_file}&
