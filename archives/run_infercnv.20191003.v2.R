# Yige Wu @WashU Oct 2019
## for running inferCNV using integrated seruat object (raw count too big)

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# input gene order file ---------------------------------------------------
tab_gene_order_ensembl = read.delim(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/copy_number/run_infercnv/20191001.v1/gencode_v21_gen_pos.ensembl_gene_id.txt", header = FALSE,stringsAsFactors = FALSE)


# set parameters ----------------------------------------------------------
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

sample_id <- "CPT0086350004"

# load seruat object to get annotation files ------------------------------------------------------
int_seurat_object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/integration/integrate_seurat_objects/20190927.v1/renal_integrated.20190927.v1.RDS")
seurat_object <- Seurat::subset(object = int_seurat_object, idents = sample_id)
saveRDS(object = seurat_object, file = paste0(dir_out, sample_id, "_from_integrated_seurat_object.RDS"))

# get barcode 2 cluster table ---------------------------------------------
anno_tab <- seurat_object@meta.data[seurat_object@meta.data$orig.ident == sample_id,]
anno_tab$barcode <- rownames(anno_tab)
anno_tab <- anno_tab %>%
  select(barcode, seurat_clusters)
write.table(x = anno_tab, file = paste0(dir_out, sample_id, "_annotation_file.", run_id, ".txt"), quote = F, sep = "\t", row.names = F, col.names = F)
rm(seurat_object)
nrow(anno_tab)


# input feature gene symbol to ensembl mapping ----------------------------
feature.names = read.delim(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/Cell_Ranger/outputs/CPT0086350004/outs/raw_feature_bc_matrix/features.tsv.gz", header = FALSE,stringsAsFactors = FALSE)

# get cells 2 process -----------------------------------------------------
sample_barcodes_int <- rownames(seurat_object@meta.data)[seurat_object@meta.data$orig.ident == sample_id]
sample_barcodes_int

# get raw read count matrix -----------------------------------------------
raw_exp_mat <- seurat_object@assays$RNA@counts[, sample_barcodes_int]
dim(raw_exp_mat)

# transform the rownames to ensembl ids -----------------------------------
raw_exp_mat <- raw_exp_mat[rownames(raw_exp_mat) %in% feature.names$V2, ]
dim(raw_exp_mat)

rownames(raw_exp_mat) <- plyr::mapvalues(rownames(raw_exp_mat), from = feature.names$V2,to = feature.names$V1)
rownames(raw_exp_mat) %>% head()

# expression data - remove gene not mapped in the gene position file
missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% tab_gene_order_ensembl$V1)]
missing_sym

## only missing ~750 genes, romove them
## only keep the barcodes in the annotation files
clean_exp_mat <- raw_exp_mat[!(rownames(raw_exp_mat) %in% missing_sym), ]
rm(raw_exp_mat)
fwrite(x = as.data.frame(clean_exp_mat), file = paste0(dir_out, "raw_matrix.", sample_id, ".tsv"), quote = F, row.names = T, col.names = T, sep = "\t")

clean_exp_mat <- as.matrix(clean_exp_mat)
rownames(clean_exp_mat) %in% tab_gene_order_ensembl$V1 %>% all # Check if only gene with position left. Should be TRUE
rownames(clean_exp_mat) %>% duplicated %>% any # Check if duplicate still exist. Should be FALSE

# get ref group names -----------------------------------------------------
all_group_names <- as.vector(unique(anno_tab$group))
all_group_names
tumor_group_names <- c("C1", "C3", "C4", "C6", "C7", "C8")
ref_group_names <- all_group_names[!(all_group_names %in% tumor_group_names)]
ref_group_names

# run infercnv ------------------------------------------------------------
infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix=clean_exp_mat,
                                    annotations_file="./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/copy_number/run_infercnv/20191002.v1/annotation_file.20191002.v1.txt",
                                    delim="\t",
                                    gene_order_file="./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/copy_number/run_infercnv/20191001.v1/gencode_v21_gen_pos.ensembl_gene_id.txt",
                                    ref_group_names= ref_group_names)

