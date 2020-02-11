#!/bin/bash

docker_image_name=sridnona/seurat_docker
dir_base=/diskmnt/Projects/ccRCC_scratch/
dir_log=${dir_base}"Logs/"
mkdir -p ${dir_log}
path_log_file=${dir_log}$0.$(date +%Y%m%d%H%M%S).log

docker run -v ${dir_base}:${dir_base} ${docker_image_name} Rscript $1 &> ${path_log_file}&
echo "Log file is at "${path_log_file}
