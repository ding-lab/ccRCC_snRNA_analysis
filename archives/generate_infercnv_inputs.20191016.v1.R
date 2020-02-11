# Yige Wu @WashU Oct 2019
## for running inferCNV using integrated seruat object (raw count too big)

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set input & output directory ----------------------------------------------------
dir_infercnv_out <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/inferCNV/outputs/"
dir.create(dir_infercnv_out)
dir_int_obj <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/integration/"

# set run id ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set integration id and output directory------------------------------------------------------
integration_id <- "20191015.v1"
dir_infercnv_out_by_int <- paste0(dir_infercnv_out, "integration.", integration_id, "/")
dir.create(dir_infercnv_out_by_int)


# input cell type assignment table ----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191016.v2.tsv", data.table = F)

# set tumor groups --------------------------------------------------------
tumor_group_names <- cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Malignant == "Yes"]
tumor_group_names <- as.character(tumor_group_names)
tumor_group_names

ref_group_names <- cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Malignant == "No"]
ref_group_names
# [1]  3  4  5  8  9 10 12 13 14 15 16

# input integrated seurat object ------------------------------------------
int_seurat_object <- readRDS("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/integration/integrate_seurat_objects/20191015.v1/Renal_Integrated.20191015.v1.RDS")
int_seurat_object@meta.data$orig.ident %>% unique()

# [1] "CPT0075140002" "CPT0001260013" "CPT0086350004" "CPT0001180011" "CPT0025890002" "CPT0010110013" "CPT0019130004"

for (snRNA_aliquot_id in unique(int_seurat_object@meta.data$orig.ident)) {
  dir_infercnv_out_by_int_sample <- paste0(dir_infercnv_out_by_int, snRNA_aliquot_id, "/")
  dir.create(dir_infercnv_out_by_int_sample)
  
  if (!file.exists(paste0(dir_infercnv_out_by_int_sample, snRNA_aliquot_id, ".Ref_Cluster_Names.", run_id, ".txt"))) {
    # subset integrated object by sample ------------------------------------------------------
    seurat_object <- Seurat::SubsetData(int_seurat_object, idents = snRNA_aliquot_id)
    
    ref_group_names2p <- intersect(ref_group_names, as.character(unique(seurat_object@meta.data$seurat_clusters)))
    
    sink(paste0(dir_infercnv_out_by_int_sample, snRNA_aliquot_id, ".Ref_Cluster_Names.", run_id, ".txt"))
    cat(paste(ref_group_names2p, collapse = ","))
    sink()
  }
}

for (snRNA_aliquot_id in unique(int_seurat_object@meta.data$orig.ident)) {
  dir_infercnv_out_by_int_sample <- paste0(dir_infercnv_out_by_int, snRNA_aliquot_id, "/")
  dir.create(dir_infercnv_out_by_int_sample)
  
  if (!file.exists(paste0(dir_infercnv_out_by_int_sample, snRNA_aliquot_id, ".RNA_Count.", run_id, ".tsv"))) {
    # subset integrated object by sample ------------------------------------------------------
    seurat_object <- Seurat::SubsetData(int_seurat_object, idents = snRNA_aliquot_id)
    
    # get barcode 2 cluster table ---------------------------------------------
    anno_tab <- seurat_object@meta.data[seurat_object@meta.data$orig.ident == snRNA_aliquot_id,]
    anno_tab$barcode <- rownames(anno_tab)
    anno_tab <- anno_tab %>%
      select(barcode, seurat_clusters)
    nrow(anno_tab)
    write.table(x = anno_tab, file = paste0(dir_infercnv_out_by_int_sample, snRNA_aliquot_id, ".Barcode_Annotation.", run_id, ".txt"), quote = F, sep = "\t", row.names = F, col.names = F)
    
    # get cells 2 process -----------------------------------------------------
    sample_barcodes_int <- anno_tab$barcode
    sample_barcodes_int
    
    # get raw read count matrix -----------------------------------------------
    raw_exp_mat <- seurat_object@assays$RNA@counts[, sample_barcodes_int]
    dim(raw_exp_mat)
    
    # expression data - remove gene not mapped in the gene position file
    tab_gene_order_symbol = read.delim(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/inferCNV/inputs/gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt", header = FALSE,stringsAsFactors = FALSE)
    missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% tab_gene_order_symbol$V1)]
    missing_sym
    
    ## only missing ~750 genes, romove them
    ## only keep the barcodes in the annotation files
    clean_exp_mat <- raw_exp_mat[!(rownames(raw_exp_mat) %in% missing_sym), ]
    rm(raw_exp_mat)
    clean_exp_mat <- as.matrix(clean_exp_mat)
    write.table(x = clean_exp_mat, file = paste0(dir_infercnv_out_by_int_sample, snRNA_aliquot_id, ".RNA_Count.", run_id, ".tsv"), quote = F, row.names = T, col.names = T, sep = "\t")
    
    # check if the matrix is good for running ---------------------------------
    rownames(clean_exp_mat) %in% tab_gene_order_symbol$V1 %>% all # Check if only gene with position left. Should be TRUE
    rownames(clean_exp_mat) %>% duplicated %>% any # Check if duplicate still exist. Should be FALSE
  }
  
}

