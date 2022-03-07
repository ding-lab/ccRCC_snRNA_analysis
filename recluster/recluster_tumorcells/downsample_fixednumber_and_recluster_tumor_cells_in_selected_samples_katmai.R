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
## specify the minical number of cells to subsample to
min_tumorcells <- 2000

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

  ## remove doublets
  ## take out the doublets
  barcode2scrublet_tmp_df <- barcode2scrublet_df %>%
    filter(Aliquot_WU == easy_id_tmp) %>%
    filter(Barcode %in% rownames(srat@meta.data)) %>%
    filter(predicted_doublet)
  barcodes_doublet <- barcode2scrublet_tmp_df$Barcode
  barcodes_keep <- colnames(srat); length(barcodes_keep)
  barcodes_keep <- barcodes_keep[!(barcodes_keep %in% barcodes_doublet)]; length(barcodes_keep)
  
  if (length(barcodes_keep) >= min_tumorcells) {
    
    ## check if the reclustered object has been saved for this aliquot
    file2write <- paste0(dir_out, easy_id_tmp, ".", min_tumorcells, "tumorcellreclustered.", run_id, ".RDS")
    if (!file.exists(file2write)) {
      barcodes_keep2 <- sample(x = barcodes_keep, size = min_tumorcells)
      
      ## subset data
      print(dim(srat))
      srat.new <- subset(srat, cells = barcodes_keep2)
      rm(srat)
      print(dim(srat.new))
      
      ## Run the standard workflow for clustering and visualization
      srat.new <- ScaleData(srat.new, features = rownames(srat.new@assays$RNA@counts))
      srat.new <- FindVariableFeatures(object = srat.new, selection.method = "vst", nfeatures = 2000)
      
      ## RunPCA
      print(num_pc)
      srat.new <- RunPCA(srat.new, npcs = num_pc, verbose = FALSE)
      srat.new <- RunUMAP(srat.new, reduction = "pca", dims = 1:num_pc)
      srat.new <- FindNeighbors(srat.new, reduction = "pca", dims = 1:num_pc, force.recalc = T)
      srat.new <- FindClusters(srat.new, resolution = 0.5)
      # srat.new <- FindClusters(srat.new, resolution = 1.0)
      saveRDS(object = srat.new, file = file2write, compress = T)
      print(paste0("Finished reclustering ", easy_id_tmp, "!"))
    } else {
      srat.new <- readRDS(file = file2write)
    }
    ## extract current meta data
    barcodes_tmp_df <- FetchData(object = srat.new, vars = c("UMAP_1", "UMAP_2", "orig.ident", "seurat_clusters"))
    barcodes_tmp_df$barcode <- rownames(barcodes_tmp_df)
    barcodes_tmp_df$easy_id <- easy_id_tmp
    barcodes_tmp_df <- barcodes_tmp_df %>%
      arrange(seurat_clusters)
    ## bind with the super table
    barcodes_umapdata_df <- rbind(barcodes_tmp_df, barcodes_umapdata_df)
    print(paste0("Finished FetchData for ", easy_id_tmp))
    
    ## plot
    p <- ggplot()
    p <- p + geom_point(data = barcodes_tmp_df, 
                        mapping = aes(x = UMAP_1, y = UMAP_2, color = factor(seurat_clusters)),
                        alpha = 1, size = 0.05)
    p <- p + scale_color_manual(values = sort(colors_cluster[unique(barcodes_tmp_df$seurat_clusters)]))
    p <- p + ggtitle(label = paste0(easy_id_tmp, " tumor cells subclusters (Seurat)"), subtitle = paste0("Down-sampling ", min_tumorcells, " cells"))
    p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = NULL, nrow = 1))
    p <- p + theme_void()
    p <- p + theme(legend.position = "bottom")
    p <- p + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15))
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank())
    ## save plot
    png2write <- paste0(dir_out, easy_id_tmp, ".png")
    png(filename = png2write, width = 900, height = 1000, res = 150)
    print(p)
    dev.off()
    print(paste0("Finished plotting for ", easy_id_tmp))
  }
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "UMAPData.", min_tumorcells, "TumorCellReclustered.", run_id, ".tsv")
write.table(x = barcodes_umapdata_df, file = file2write, sep = '\t', quote = F, row.names = F)
print("Finished write.table!")
