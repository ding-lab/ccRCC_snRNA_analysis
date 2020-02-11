# Yige Wu @WashU Feb 2020
## for running inferCNV using individual seruat objects

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set input & output directory ----------------------------------------------------
dir_infercnv <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/"
dir_infercnv_inputs <- paste0(dir_infercnv, "inputs/")
dir.create(dir_infercnv_inputs)
dir_infercnv_counts <- paste0(dir_infercnv_inputs, "raw_counts_matrix/")
dir.create(dir_infercnv_counts)
dir_infercnv_annotation <- paste0(dir_infercnv_inputs, "annotations_file/")
dir.create(dir_infercnv_annotation)

# set run id ----------------------------------------------------------
version_tmp <- 1
# run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
run_id <- "20200207.v1"

# set output directory------------------------------------------------------
dir_infercnv_annotation_out <- paste0(dir_infercnv_annotation, "Individual.", run_id, "/")
dir.create(dir_infercnv_annotation_out)

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200207.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(FACS == "") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object

# write annotation files -------------------------------------------------
for (snRNA_aliquot_id_tmp in unique(seurat_summary2process$Aliquot)) {
  path_annotations_file_out <- paste0(dir_infercnv_annotation_out, snRNA_aliquot_id_tmp, ".Barcode_Annotation.txt")
  
  if (!file.exists(path_annotations_file_out)) {
    ## input the seurat object
    seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
    seurat_object <- readRDS(file = seurat_obj_path)
    
    ## get barcode to cluster annotation from the meta data
    anno_tab <- seurat_object@meta.data
    anno_tab$barcode <- rownames(anno_tab)
    anno_tab <- anno_tab %>%
      select(barcode, seurat_clusters)
    nrow(anno_tab)
    write.table(x = anno_tab, file = path_annotations_file_out, quote = F, sep = "\t", row.names = F, col.names = F)
  }
}