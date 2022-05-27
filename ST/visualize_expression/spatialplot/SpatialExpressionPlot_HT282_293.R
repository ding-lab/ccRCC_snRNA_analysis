## THE START FROM SCARTCH VERSION
library(tidyverse)
library(Seurat)
library(patchwork)

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2"
)
## if package is not installed, it will be installed here
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# Load
ffpe_path = '~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/Spatial_transcriptomics/output_FFPE'
stobjs = list()
stobjs$ht282 = Load10X_Spatial(str_glue('{ffpe_path}/HT282N1/H3/HT282N1-S1H3Fs4U1Bp1/outs'))
stobjs$ht293 = Load10X_Spatial(str_glue('{ffpe_path}/HT293N1/H3/HT293N1-S1H3Fs1U1Bp1/outs'))

#######################################
### Quick QC [Optional]
## Note - FFPE no MT genes. MT zero
stobjs = map(stobjs, function(obj){
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  return(obj)
  })

# # Visualize QC metrics
# ## 1. UMI Count
# p_ncount = map(stobjs, function(obj){
#   plot1 <- VlnPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3, pt.size = 0.1) & NoLegend()
#   plot2 <- SpatialPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3) & theme(legend.position = "right")
#   wrap_plots(plot1, plot2, ncol=1)
# }) %>% setNames(names(stobjs))
# # Check 
# p_ncount$ht282
# p_ncount$ht293

#######################################
# SCTransform
# https://satijalab.org/seurat/articles/spatial_vignette.html
stobjs = map(stobjs, SCTransform, assay = "Spatial", verbose = T)

#######################################
# Spatial Plot
genes_check = c("LAMA1", "LAMA3", "COL5A1", "COL8A1", "COL28A1", "CD70", "CD27", "OSM", "OSMR", "IL6ST")
genes_check = c("COL5A2", "COL6A3", "COL16A1", "COL8A1", "COL5A2", "COL28A1", "COL4A1", "COL7A1", "COL6A2", "COL4A2", "COL18A1", "COL17A1", "COL5A1",
                "COL12A1", "COL1A1", "COL4A3", "COL27A1", "COL6A1", "COL10A1", "COL4A4", "COL4A4",
                "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4")

file2write <- paste0(dir_out, "HT283.collagens.", "png")
png(file2write, width = 2000, height = 3000, res = 150)
SpatialPlot(stobjs$ht282, features = genes_check, stroke= NA, image.alpha=0, ncol = 5)
dev.off()

file2write <- paste0(dir_out, "HT293.collagens.", "png")
png(file2write, width = 2000, height = 3000, res = 150)
SpatialPlot(stobjs$ht293, features = genes_check, stroke= NA, image.alpha=0, ncol = 5)
dev.off()

for (gene_plot in c("CP", "COL4A1", "COL4A3", "LOX")) {
  file2write <- paste0(dir_out, "HT283.", gene_plot, ".png")
  png(file2write, width = 800, height = 800, res = 150)
  # print(SpatialPlot(stobjs$ht282, features = gene_plot, stroke= NA, image.alpha=0, ncol = 5))
  SpatialPlot(stobjs$ht282, features = gene_plot, stroke= NA, image.alpha=0, ncol = 5)  
  dev.off()
  
  # file2write <- paste0(dir_out, "HT293.collagens.", "png")
  # png(file2write, width = 2000, height = 3000, res = 150)
  # SpatialPlot(stobjs$ht293, features = genes_check, stroke= NA, image.alpha=0, ncol = 5)
  # dev.off()
}
