# Yige Wu @WashU Feb 2020
## for integrating tumor and normal snRNA datasets belong to the same patient

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# set case ids to be processed --------------------------------------------
case_id <- c("C3N-01200")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200207.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Case %in% case_id) %>%
  filter(FACS == "") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object

# get aliquot ids to be processed -----------------------------------------
snRNA_aliquot_ids <- seurat_summary2process$Aliquot
snRNA_aliquot_ids

# Input seurat objects into a list -----------------------------------------------------
renal.list <- list()
for (i in 1:nrow(seurat_summary2process)) {
  sample_id_tmp <- seurat_summary2process$Aliquot[i]
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[i]
  seurat_obj <- readRDS(file = seurat_obj_path)
  renal.list[[sample_id_tmp]] <- seurat_obj
  
}
length(renal.list)
rm(seurat_obj)


# Merge Data -----------------------------------------------------------
renal.merged <- merge(x = renal.list[[1]], y = renal.list[2:length(renal.list)], project = "Integrated_without_anchor")
# Warning message:
#   In CheckDuplicateCellNames(object.list = objects) :
#   Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.

# Run the standard workflow for visualization and clustering ------------
# The variable features of this assay are automatically set during IntegrateData
# renal.merged <- ScaleData(renal.merged, verbose = F) 
renal.merged <- SCTransform(renal.merged, vars.to.regress = c("nCount_RNA","percent.mito"), return.only.var.genes = F)
renal.merged <- RunPCA(renal.merged, npcs = 30, verbose = FALSE)
renal.merged <- RunUMAP(renal.merged, reduction = "pca", dims = 1:30)

renal.merged <- FindNeighbors(renal.merged, reduction = "pca", dims = 1:30)
renal.merged <- FindClusters(renal.merged, resolution = 0.5)

saveRDS(object = renal.merged, file = paste0(dir_out, case_id, ".Tumor_Normal.Integration_without_anchor.", run_id, ".RDS"))
