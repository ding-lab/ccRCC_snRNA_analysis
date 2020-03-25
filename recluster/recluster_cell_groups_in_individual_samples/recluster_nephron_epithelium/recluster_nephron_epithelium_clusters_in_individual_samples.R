# Yige Wu @WashU Feb 2020
## for each individual sample, isolating cells within the immune cell clusters and re-do clustering

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
## input seurat object master list
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200207.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object
## input the cluster to cell type table
cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20200207.v2.tsv", data.table = F)

# run reclustering by each aliquot ----------------------------------------
for (snRNA_aliquot_id_tmp in seurat_summary2process$Aliquot) {
  ## create output directory by aliquot
  dir_out1 <- paste0(dir_out, snRNA_aliquot_id_tmp, "/")
  dir.create(dir_out1)
  
  ## check if the reclustered object has been saved for this aliquot
  file2write <- paste0(dir_out1, snRNA_aliquot_id_tmp, ".nephron_epithelium_reclustered.", run_id, ".RDS")
  if (!file.exists(file2write)) {
    ## input individually processed seurat object
    seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
    seurat_obj_path
    seurat_object <- readRDS(file = seurat_obj_path)
    
    ## get the immune clusters for this aliquot
    clusters2process <- cluster2celltype_df$Cluster[cluster2celltype_df$Aliquot == snRNA_aliquot_id_tmp & cluster2celltype_df$Most_Enriched_Cell_Group == "Nephron_Epithelium"]
    clusters2process
    
    ## subset data
    new.object <- subset(seurat_object, idents = clusters2process)
    rm(seurat_object)
    
    ## Run the standard workflow for clustering and visualization
    new.object <- FindVariableFeatures(object = new.object, selection.method = "vst", nfeatures = 2000)
    new.object <- ScaleData(new.object, features = rownames(new.object@assays$RNA@counts)) 
    new.object <- RunPCA(new.object, npcs = 50, verbose = FALSE)
    new.object <- RunUMAP(new.object, reduction = "pca", dims = 1:30)
    new.object <- FindNeighbors(new.object, reduction = "pca", dims = 1:30, force.recalc = T)
    new.object <- FindClusters(new.object, resolution = 0.5)
    saveRDS(object = new.object, file = file2write)
  }
}
