# Yige Wu @WashU Dec 2020
## for running inferCNV using individual seruat objects

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# set parameters ----------------------------------------------------------
## set output directory
dir_infercnv <- "./Resources/snRNA_Processed_Data/InferCNV/"
dir_infercnv_inputs <- paste0(dir_infercnv, "inputs/")
dir.create(dir_infercnv_inputs)
dir_infercnv_annotation <- paste0(dir_infercnv_inputs, "annotations_file/")
dir.create(dir_infercnv_annotation)
dir_infercnv_annotation_out <- paste0(dir_infercnv_annotation, "Individual.", run_id, "/")
dir.create(dir_infercnv_annotation_out)

# input barcode 2 cell type table -----------------------------------------
## input barcode to cell type info
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20201130.v1/31Aliquot.Barcode2CellType.20201130.v1.tsv", data.table = F)

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
