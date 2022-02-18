# Yige Wu @WashU Nov 2020
## these samples are newly added samples so tumor-cell-only clustering take place after taking out the doublets
## the previous tumor-cell-only clustering did not take out the doublets for minimal manual work (assignment of the tumor clusters)

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
dir_tmp <- getwd()
print(dir_tmp)
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
srat_paths <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20210729.v1/Seurat_Object_Paths.20210729.v1.tsv")
### filter out normal sample
srat_paths <- srat_paths %>%
  filter(Sample_Type == "Tumor")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input the cell to cell type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
## input doublet prediction
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)
### specify aliquots to process
# aliquots2process <- unique(srat_paths$Aliquot)
aliquots2process <- c("CPT0012550012")
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
    print(dim(srat))
    
    ## get the tumor cell barcodes for this aliquot
    barcodes2process <- barcode2celltype_df$individual_barcode[barcode2celltype_df$orig.ident == aliquot_tmp & barcode2celltype_df$Cell_group5 == "Tumor cells"]
    
    ## take out the doublets
    barcode2scrublet_tmp_df <- barcode2scrublet_df %>%
      filter(Aliquot == aliquot_tmp) %>%
      filter(Barcode %in% rownames(srat@meta.data)) %>%
      filter(!predicted_doublet) %>%
      filter(Barcode %in% barcodes2process)
    barcodes_keep <- barcode2scrublet_tmp_df$Barcode
    
    ## subset data
    srat.new <- subset(srat, cells = barcodes_keep)
    rm(srat)
    print(dim(srat.new))
    
    ## Run the standard workflow for clustering and visualization
    srat.new <- ScaleData(srat.new, features = rownames(srat.new@assays$RNA@counts))
    srat.new <- FindVariableFeatures(object = srat.new, selection.method = "vst", nfeatures = 2000)
    
    ## RunPCA
    srat.new <- RunPCA(srat.new, npcs = num_pc, verbose = FALSE)
    srat.new <- RunHarmony(srat.new)
    
    srat.new <- RunUMAP(srat.new, reduction = "harmony")
    srat.new <- FindNeighbors(srat.new, reduction = "harmony", dims = 1:num_pc, force.recalc = T)
    srat.new <- FindClusters(srat.new, resolution = 0.5)
    # srat.new <- FindClusters(srat.new, resolution = 1.0)
    saveRDS(object = srat.new, file = file2write, compress = T)
    
    p <- DimPlot(object = srat.new, reduction = "harmony", pt.size = .1, do.return = TRUE)
    file_plot<- paste(dir_out, easyid, ".dimplot.png")
    png(file_plot, width = 600, height = 600, res = 150)
    print(p)
    dev.off()
    
    print(paste0("Finished ", aliquot_tmp, "!"))
  }
}

print("Finished all!")
