#!/bin/bash

docker_image_name=sridnona/seurat_docker
dir_base=/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/
#dir_log=${dir_base}"Logs/"
#mkdir -p ${dir_log}
#path_log_file=${dir_log}$0.$(date +%Y%m%d%H%M%S).log

## make timestamp
timestamp=$(date +%Y%m%d%H%M%S)

## make sure the 2nd argument is the output directory
mkdir -p $2

## make directory for the run info
dir_run_info=$2"/"${timestamp}"/"
mkdir -p ${dir_run_info}

## make log file
path_log_file=$2"/"${timestamp}".log"
path_log_file=${dir_run_info}${timestamp}".log"

## save the code for running docker
path_main_file=$2"/"${timestamp}".main.sh"
path_main_file=${dir_run_info}${timestamp}".main.sh"
echo "docker run -v ${dir_base}:${dir_base} ${docker_image_name} Rscript $1 ${@:2} &> ${path_log_file}&" > ${path_main_file}

## make a copy of the rscript
cp $1 ${dir_run_info}

## run docker
docker run -v ${dir_base}:${dir_base} ${docker_image_name} Rscript $1 ${@:2} &> ${path_log_file}&

## print path to the current run info
echo ""
echo "Output directory is at "$2
echo "Code for the docker run is at "${path_main_file}
echo "Log file is at "${path_log_file}
echo ""
