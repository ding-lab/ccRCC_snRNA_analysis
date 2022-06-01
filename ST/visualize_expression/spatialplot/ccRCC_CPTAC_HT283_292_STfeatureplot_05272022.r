## Spatial Analysis after Label transfer 
## Spatial Plot HT282 and 293
## R version
library('patchwork')
library('tidyverse')
library("Seurat")

setwd('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/RCC')

today_date = format(Sys.time(), '%m%d%Y')
ptm = Sys.time()

####################################### 
## Load and process
source('script/ccRCC_loadST_08092021_wo_RCTD.R')

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
p1 = SpatialPlot(st_htan, features = genes_check, stroke = NA, image.alpha = 0)

##############
# v2: Plot with shared feature scale
## Plot function
SpatialPlotSharedScale = function(obj, feature, range_ratio = 0.9, ...){
    exp_max = FetchData(obj, feature) %>% max() %>% "*"(range_ratio)
    SpatialPlot(st_htan, features = feature, ...) & scale_fill_gradientn(colors = Seurat:::SpatialColors(100), limit = c(0, exp_max))
}
SpatialPlotSharedScaleMultipleGenes = function(obj, features, range_ratio = 0.9, ...){
    p_list = map(features, function(gene) SpatialPlotSharedScale(obj, feature = gene, range_ratio=range_ratio, ...))
    wrap_plots(p_list, ncol =1 )
}

p2 = SpatialPlotSharedScaleMultipleGenes(st_htan, features = genes_check, stroke = NA, image.alpha = 0)

##############               
## Save
dir.create('figure/STExp/05272022/')
pdf('figure/STExp/05272022/ccRCC_HT282_293_multiple_genes.pdf', width = 8, height = 25)
print(p1)
print(p2)
dev.off()