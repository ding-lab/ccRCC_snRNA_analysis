## Spatial Analysis after Label transfer 
## Spatial Plot HT282 and 293
## R version
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
# dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "patchwork",
  "tidyverse"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

####################################### 
## Load and process
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/RCC/script/ccRCC_loadST_08092021_wo_RCTD.R')

## Subset HT282 and HT293
ht_ids = st_merge@meta.data %>% filter(str_detect(orig.ident, 'HT282|HT293')) %>% rownames
st_htan = subset(st_merge, cells = ht_ids)

## Clean Image slot
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Shared_resource/script/function_seurat_janitor.R')
st_htan = FixSeuratImage(st_htan) 


#######################################
# Plotting
# v1: Direct Plot with different scale
genes_check = c("LAMA1", "LAMA3", "COL5A1", "COL8A1", "COL28A1", "CD70", "CD27", "OSM", "OSMR", "IL6ST")
genes_check = c("CA9", "CP", "COL4A1", "OSMR", "TGM2")
genes_check = c("CA9", "CP", "COL4A1", "OSM", "TGM2")
p1 = SpatialPlot(st_htan, features = genes_check, stroke = NA, image.alpha = 0, ncol = 5)

##############
# v2: Plot with shared feature scale
## Plot function
SpatialPlotSharedScale = function(obj, feature, range_ratio = 0.9, ...){
    exp_max = FetchData(obj, feature) %>% max() %>% "*"(range_ratio)
    SpatialPlot(st_htan, features = feature, ...) & scale_fill_gradientn(colors = Seurat:::SpatialColors(100), limit = c(0, exp_max))
}
SpatialPlotSharedScaleMultipleGenes = function(obj, features, range_ratio = 0.9, ...){
    p_list = map(features, function(gene) SpatialPlotSharedScale(obj, feature = gene, range_ratio=range_ratio, ...))
    wrap_plots(p_list, nrow = 5)
}

p2 = SpatialPlotSharedScaleMultipleGenes(st_htan, features = genes_check, stroke = NA, image.alpha = 0)

##############               
## set output directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
source("./ccRCC_snRNA_analysis/functions.R")
dir_out_parent <- makeOutDir_katmai(path_this_script)
dir_out <- paste0(dir_out_parent, run_id, "/")
dir.create(dir_out)
## Save
pdf(paste0(dir_out, "ccRCC_HT282_293.", "different_scales.", "pdf"), width = 6, height = 15, useDingbats = F)
print(p1)
dev.off()

pdf(paste0(dir_out, "ccRCC_HT282_293.", "same_scales.", "pdf"), width = 6, height = 15, useDingbats = F)
print(p2)
dev.off()
