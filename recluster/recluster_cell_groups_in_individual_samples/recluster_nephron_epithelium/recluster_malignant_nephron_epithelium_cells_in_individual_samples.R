# Yige Wu @WashU Feb 2020
## for each individual sample, isolating tumor cells assigned from the integrated data and re-do clustering

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input seurat object paths
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_cluster/write_path_to_srats_on_box/20200219.v1/srat_Paths.20200219.v1.tsv", data.table = F)
### filter out normal sample
srat_paths <- srat_paths %>%
  filter(Sample_Type == "Tumor")
## input the cell to cell type table
barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/map_celltype_to_barcode/20200224.v1/30_aliquot_integration.barcode2celltype.20200224.v1.tsv", data.table = F)

# run reclustering by each aliquot ----------------------------------------
for (aliquot_tmp in unique(srat_paths$Aliquot)) {
  ## check if the reclustered object has been saved for this aliquot
  file2write <- paste0(dir_out, aliquot_tmp, ".malignant_nephron_epithelium.", run_id, ".RDS")
  if (!file.exists(file2write)) {
    ## input individually processed seurat object
    seurat_obj_path <- srat_paths$Path_srat[srat_paths$Aliquot == aliquot_tmp]
    seurat_obj_path
    srat <- readRDS(file = seurat_obj_path)
    
    ## get the tumor cell barcodes for this aliquot
    barcodes2process <- barcode2celltype_df$individual_barcode[barcode2celltype_df$orig.ident == aliquot_tmp & barcode2celltype_df$Most_Enriched_Cell_Group == "Nephron_Epithelium" & barcode2celltype_df$Is_Normal_Nephron_Epithelium == F]
    
    ## subset data
    srat.new <- subset(srat, cells = barcodes2process)
    rm(srat)
    
    ## Run the standard workflow for clustering and visualization
    srat.new <- FindVariableFeatures(object = srat.new, selection.method = "vst", nfeatures = 2000)
    srat.new <- ScaleData(srat.new, features = rownames(srat.new@assays$RNA@counts)) 
    srat.new <- RunPCA(srat.new, npcs = 30, verbose = FALSE)
    srat.new <- RunUMAP(srat.new, reduction = "pca", dims = 1:30)
    srat.new <- FindNeighbors(srat.new, reduction = "pca", dims = 1:30, force.recalc = T)
    srat.new <- FindClusters(srat.new, resolution = 0.5)
    saveRDS(object = srat.new, file = file2write, compress = T)
  }
}

# make table for paths to these objects -----------------------------------
new_srat_paths <- data.frame(Aliquot = unique(srat_paths$Aliquot),
                             Path_srat = paste0(dir_out, unique(srat_paths$Aliquot), ".malignant_nephron_epithelium.", run_id, ".RDS"))
write.table(x = new_srat_paths, file = paste0(dir_out, "srat_Paths", ".", "Malignant_Nephron_Epithelium", run_id, ".", "tsv"), 
            quote = F, sep = "\t", row.names = F)


