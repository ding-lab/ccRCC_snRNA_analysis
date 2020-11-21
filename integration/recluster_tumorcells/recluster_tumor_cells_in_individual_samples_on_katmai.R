# Yige Wu @WashU Nov 2020
## for each individual sample, isolating tumor cells assigned from the integrated data and re-do clustering

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
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
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

# input dependencies ------------------------------------------------------
## input seurat object paths
srat_paths <- fread(input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20200717.v1/Seurat_Object_Paths.20200717.v1.tsv", data.table = F)
### filter out normal sample
srat_paths <- srat_paths %>%
  filter(Sample_Type == "Tumor")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input the cell to cell type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20201121.v1/31Aliquot.Barcode2CellType.20201121.v1.tsv", data.table = F)
### specify aliquots to process
aliquots2process <- unique(srat_paths$Aliquot)
aliquots2process

# run reclustering by each aliquot ----------------------------------------
for (aliquot_tmp in aliquots2process) {
  easyid <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquot_tmp]
  ## check if the reclustered object has been saved for this aliquot
  file2write <- paste0(dir_out, easyid, ".tumorcellreclustered.", run_id, ".RDS")
  if (!file.exists(file2write)) {
    ## input individually processed seurat object
    seurat_obj_path <- srat_paths$Path_katmai_seurat_object[srat_paths$Aliquot == aliquot_tmp]
    seurat_obj_path
    srat <- readRDS(file = seurat_obj_path)
    
    ## get the tumor cell barcodes for this aliquot
    barcodes2process <- barcode2celltype_df$individual_barcode[barcode2celltype_df$orig.ident == aliquot_tmp & barcode2celltype_df$Cell_group5 == "Tumor cells"]
    
    ## subset data
    srat.new <- subset(srat, cells = barcodes2process)
    rm(srat)
    
    ## Run the standard workflow for clustering and visualization
    srat.new <- ScaleData(srat.new, features = rownames(srat.new@assays$RNA@counts))
    srat.new <- FindVariableFeatures(object = srat.new, selection.method = "vst", nfeatures = 2000)
    
    ## RunPCA
    srat.new <- RunPCA(srat.new, npcs = num_pc, verbose = FALSE)
    srat.new <- RunUMAP(srat.new, reduction = "pca", dims = 1:num_pc)
    srat.new <- FindNeighbors(srat.new, reduction = "pca", dims = 1:num_pc, force.recalc = T)
    srat.new <- FindClusters(srat.new, resolution = 0.5)
    # srat.new <- FindClusters(srat.new, resolution = 1.0)
    saveRDS(object = srat.new, file = file2write, compress = T)
    print(paste0("Finished all!"))
  }
}

print("Finished all!")
