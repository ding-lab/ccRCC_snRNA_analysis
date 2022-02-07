# Yige Wu @WashU Dec 2020
## for running inferCNV using individual seruat objects

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
dir_infercnv_annotation_out <- paste0(dir_infercnv_annotation, "run.", "20220201", "/")
dir.create(dir_infercnv_annotation_out)

# input barcode 2 cell type table -----------------------------------------
## input barcode to cell type info
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
## input doublet info
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# write annotation files -------------------------------------------------
aliquots_nat <- metadata_df$Aliquot.snRNA[metadata_df$Sample_Type == "Normal"]
for (snRNA_aliquot_id_tmp in unique(barcode2celltype_df$orig.ident)) {
  path_annotations_file_out <- paste0(dir_infercnv_annotation_out, snRNA_aliquot_id_tmp, ".Barcode_Annotation.txt")
  
  if (!file.exists(path_annotations_file_out)) {
    ## get barcode to cluster annotation from the meta data
    if (!(snRNA_aliquot_id_tmp %in% aliquots_nat)) {
      anno_tab <- barcode2celltype_df %>%
        filter(orig.ident == snRNA_aliquot_id_tmp) %>%
        mutate(infercnv_group = ifelse(Cell_group5 %in% c("Tumor cells", "Unknown"), "Obs", "Ref")) %>%
        select(individual_barcode, infercnv_group)
    } else {
      anno_tab <- barcode2celltype_df %>%
        filter(orig.ident == snRNA_aliquot_id_tmp) %>%
        mutate(infercnv_group = ifelse(Cell_group5 %in% c("Normal epithelial cells", "Tumor cells", "Unknown"), "Obs", "Ref")) %>%
        select(individual_barcode, infercnv_group)
    }

    
    if (snRNA_aliquot_id_tmp %in% barcode2scrublet_df$Aliquot) {
      barcodes_doublet <- barcode2scrublet_df$Barcode[barcode2scrublet_df$Aliquot == snRNA_aliquot_id_tmp & barcode2scrublet_df$predicted_doublet]
      anno_tab <- anno_tab %>%
        filter(!(individual_barcode %in% barcodes_doublet))
    }
    
    write.table(x = anno_tab, file = path_annotations_file_out, quote = F, sep = "\t", row.names = F, col.names = F)
  }
}
