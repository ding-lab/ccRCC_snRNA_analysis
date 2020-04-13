# Yige Wu @WashU Apr 2020
## for running inferCNV using individual seruat objects

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
run_id <- "20200413.v1"
## set output directory
dir_infercnv <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/"
dir_infercnv_inputs <- paste0(dir_infercnv, "inputs/")
dir.create(dir_infercnv_inputs)
dir_infercnv_annotation <- paste0(dir_infercnv_inputs, "annotations_file/")
dir.create(dir_infercnv_annotation)
dir_infercnv_annotation_out <- paste0(dir_infercnv_annotation, "Individual.", run_id, "/")
dir.create(dir_infercnv_annotation_out)

# input barcode 2 cell type table -----------------------------------------
## input barcode to cell type info
barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200406.v1/30_aliquot_integration.barcode2celltype.20200406.v1.tsv", data.table = F)

# write annotation files -------------------------------------------------
for (snRNA_aliquot_id_tmp in unique(barcode2celltype_df$orig.ident)) {
  path_annotations_file_out <- paste0(dir_infercnv_annotation_out, snRNA_aliquot_id_tmp, ".Barcode_Annotation.txt")
  
  if (!file.exists(path_annotations_file_out)) {
    ## get barcode to cluster annotation from the meta data
    anno_tab <- barcode2celltype_df %>%
      filter(orig.ident == snRNA_aliquot_id_tmp) %>%
      mutate(infercnv_group = Cell_type.shorter) %>%
      select(individual_barcode, infercnv_group)
    
    anno_tab$infercnv_group <- gsub(x = anno_tab$infercnv_group, pattern = '[/ +]', replacement = ".")
    anno_tab$infercnv_group <- gsub(x = anno_tab$infercnv_group, pattern = '\\-', replacement = ".")
    ## group ref cell types with too few cells into one group
    anno_count_df <- data.frame(table(anno_tab$infercnv_group))
    anno_count_regroup_df <- anno_count_df %>%
      filter(Freq < 10)
    anno_tab$infercnv_group[anno_tab$infercnv_group %in% anno_count_regroup_df$Var1] <- "Ref_SmallNumber"
    
    nrow(anno_tab)
    write.table(x = anno_tab, file = path_annotations_file_out, quote = F, sep = "\t", row.names = F, col.names = F)
  }
}
