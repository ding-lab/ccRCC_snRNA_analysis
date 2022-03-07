# Yige Wu @WashU Feb 2022

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
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
library(ggplot2)

# input dependencies ------------------------------------------------------
## input seurat object paths
srat_paths_df <- fread(data.table = F, input = "./Data_Freezes/V2/snRNA/Tumor_Cell_Reclustered/Paths_TumorCellOnlyReclustered_SeuratObject.20210805.v1.tsv")
## input the scrublet output
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)
## specify the minical number of cells to downsample to
num_pc <- 50
# run reclustering by each aliquot ----------------------------------------
barcodes_umapdata_df <- NULL
pct_df <- NULL

### make colors
colors_cluster <- Polychrome::dark.colors(n = 24)
names(colors_cluster) <- 0:23

for (easy_id_tmp in srat_paths_df$Aliquot.snRNA.WU) {
  ## input individually processed seurat object
  seurat_obj_path <- srat_paths_df$Path_katmai[srat_paths_df$Aliquot.snRNA.WU == easy_id_tmp]
  srat <- readRDS(file = seurat_obj_path)
  number_tumorcells <- ncol(srat)
  ## check if the reclustered object has been saved for this aliquot
  file2write <- paste0(dir_out, easy_id_tmp, ".PC", num_pc, "tumorcellreclustered.", run_id, ".RDS")
  if (!file.exists(file2write)) {
    ## remove doublets
    ## take out the doublets
    barcode2scrublet_tmp_df <- barcode2scrublet_df %>%
      filter(Aliquot_WU == easy_id_tmp) %>%
      filter(Barcode %in% rownames(srat@meta.data)) %>%
      filter(!predicted_doublet)
    barcodes_keep <- barcode2scrublet_tmp_df$Barcode
    
    ## subset data
    print(dim(srat))
    srat <- subset(srat, cells = barcodes_keep)
    print(dim(srat))
    
    ## RunPCA
    srat <- RunPCA(srat, npcs = num_pc, verbose = FALSE)
    srat <- RunUMAP(srat, reduction = "pca", dims = 1:num_pc)
    srat <- FindNeighbors(srat, reduction = "pca", dims = 1:num_pc, force.recalc = T)
    srat <- FindClusters(srat, resolution = 0.5)
    # srat.new <- FindClusters(srat.new, resolution = 1.0)
    saveRDS(object = srat, file = file2write, compress = T)
    print(paste0("Finished reclustering ", easy_id_tmp, "!"))
  } else {
    srat <- readRDS(file = file2write)
  }
  # Determine percent of variation associated with each PC
  pct <- srat[["pca"]]@stdev / sum(srat[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  pct_tmp_df <- data.frame(easy_id = rep(easy_id_tmp, length(pct)), rank_pc = 1:length(pct), pct = pct, cumu_pct = cumu)
  pct_df <- rbind(pct_tmp_df, pct_df)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  print(co1)
  
  ## extract current meta data
  barcodes_tmp_df <- FetchData(object = srat, vars = c("UMAP_1", "UMAP_2", "orig.ident", "seurat_clusters"))
  barcodes_tmp_df$barcode <- rownames(barcodes_tmp_df)
  barcodes_tmp_df$easy_id <- easy_id_tmp
  
  ## bind with the super table
  barcodes_umapdata_df <- rbind(barcodes_tmp_df, barcodes_umapdata_df)
  print(paste0("Finished FetchData for ", easy_id_tmp))
  
  ## plot
  p <- ggplot()
  p <- p + geom_point(data = barcodes_tmp_df, 
                      mapping = aes(x = UMAP_1, y = UMAP_2, color = factor(seurat_clusters)),
                      alpha = 1, size = 0.05)
  p <- p + scale_color_manual(values = colors_cluster[sort(unique(barcodes_tmp_df$seurat_clusters))])
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p <- p + ggtitle(label = paste0(" Tumor-cells subclusters for sample ", easy_id_tmp), subtitle = paste0("Using 50 PCs"))
  p <- p + theme(legend.position = "bottom")
  p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = NULL))
  p
  file2write <- paste0(dir_out, easy_id_tmp, ".png")
  png(filename = file2write, width = 600, height = 700, res = 150)
  print(p)
  dev.off()
  print(paste0("Finished plotting for ", easy_id_tmp))
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "UMAPData.PC", num_pc, "TumorCellReclustered.", run_id, ".tsv")
write.table(x = barcodes_umapdata_df, file = file2write, sep = '\t', quote = F, row.names = F)
print("Finished write.table for barcodes_umapdata_df!")

file2write <- paste0(dir_out, "Pct.SD.byPC.", run_id, ".tsv")
write.table(x = pct_df, file = file2write, sep = '\t', quote = F, row.names = F)
print("Finished write.table for pct_df!")
