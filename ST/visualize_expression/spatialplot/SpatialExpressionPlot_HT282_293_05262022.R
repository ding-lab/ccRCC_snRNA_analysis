## THE START FROM SCARTCH VERSION
library(tidyverse)
library(Seurat)
library(patchwork)

# Local 
setwd('~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/')

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

# Visualize QC metrics
## 1. UMI Count
p_ncount = map(stobjs, function(obj){
  plot1 <- VlnPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3, pt.size = 0.1) & NoLegend()
  plot2 <- SpatialPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3) & theme(legend.position = "right")
  wrap_plots(plot1, plot2, ncol=1)
}) %>% setNames(names(stobjs))
# Check 
p_ncount$ht282
p_ncount$ht293

#######################################
# SCTransform
# https://satijalab.org/seurat/articles/spatial_vignette.html
stobjs = map(stobjs, SCTransform, assay = "Spatial", verbose = T)

#######################################
# Spatial Plot
genes_check = c("LAMA1", "LAMA3", "COL5A1", "COL8A1", "COL28A1", "CD70", "CD27", "OSM", "OSMR", "IL6ST")
genes_check = c("COL5A2", "COL6A3", "COL16A1", "COL8A1", "COL5A2", "COL28A1", "COL4A1", "COL7A1", "COL6A2", "COL4A2", "COL18A1", "COL17A1", "COL5A1",
                "COL12A1", "COL1A1", "COL4A3")

SpatialPlot(stobjs$ht282, features = genes_check, stroke= NA, image.alpha=0)
SpatialPlot(stobjs$ht293, features = genes_check, stroke= NA, image.alpha=0)

