# Yige Wu @WashU Jun 2021

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
# dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir_out <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Data_Freezes/V2/snRNA/Individual_Seurat_Objects/"
dir.create(dir_out)
options(future.globals.maxSize = 1000 * 1024^2)

# input dependencies ------------------------------------------------------
## input seurat paths
paths_srat <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20210729.v1/Seurat_Object_Paths.20210729.v1.tsv")
## input doublet prediction
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")
## input idemta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# input the seurat object and store in a list--------------------------------------------------------
paths_rds <- NULL
easy_ids <- NULL
for (i in 1:nrow(paths_srat)) {
  sample_id_tmp <- paths_srat$Aliquot[i]
  easy_id_tmp <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == sample_id_tmp]
  path_rds_tmp <- paste0(dir_out, easy_id_tmp, ".DoubletRemoved.RDS")
  print(easy_id_tmp)
  
  if (!file.exists(path_rds_tmp)) {
    seurat_obj_path <- paths_srat$Path_katmai_seurat_object[i]
    
    seurat_obj <- readRDS(file = seurat_obj_path)
    ## take out the doublets
    barcode2scrublet_tmp_df <- barcode2scrublet_df %>%
      filter(Aliquot == sample_id_tmp) %>%
      filter(Barcode %in% rownames(seurat_obj@meta.data)) %>%
      filter(!predicted_doublet)
    barcodes_keep <- barcode2scrublet_tmp_df$Barcode
    seurat_sub_obj <- subset(x = seurat_obj, cells = barcodes_keep)
    print(dim(seurat_sub_obj))
    print("Subsetted")
    
    ## add cell type as identity
    barcode2celltype_tmp_df <- barcode2celltype_df %>%
      filter(orig.ident == sample_id_tmp)
    seurat_sub_obj@meta.data$Cell_group_w_epithelialcelltypes <- mapvalues(x = rownames(seurat_sub_obj@meta.data), from = barcode2celltype_tmp_df$individual_barcode, to = as.vector(barcode2celltype_tmp_df$Cell_group_w_epithelialcelltypes))
    Idents(seurat_sub_obj) <- "Cell_group_w_epithelialcelltypes"
    print("Added cell type")
    
    ## write output
    saveRDS(object = seurat_sub_obj, file = path_rds_tmp, compress = T)
    print("saveRDS")
  }
  paths_rds <- c(paths_rds, path_rds_tmp)
  easy_ids <- c(easy_ids, easy_id_tmp)
}
rds_list_df <- data.frame(easy_ids, paths_rds)
file2write <- paste0(dir_out, "Seurat_Object_Paths.DoubletRemoved.", run_id, ".tsv")
write.table(x = rds_list_df, file = file2write, quote = F, sep = "\t", row.names = F, col.names = F)
print("Finished All!")


