# Yige Wu @WashU Oct 2019
## for running inferCNV using integrated seruat object (raw count too big)

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# input gene order file ---------------------------------------------------
tab_gene_order_ensembl = read.delim(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/copy_number/run_infercnv/20191001.v1/gencode_v21_gen_pos.ensembl_gene_id.txt", header = FALSE,stringsAsFactors = FALSE)

# set parameters ----------------------------------------------------------
sample_ids <- c("CPT0075130004_notFACS", "CPT0086820004_notFACS", "CPT0075140002", "CPT0001260013", "CPT0086350004")

version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

sample_id <- "CPT0086350004"

# load seruat object to get annotation files ------------------------------------------------------
seurat_object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/integration/integrate_seurat_objects/20190927.v1/renal_integrated.20190927.v1.RDS")

# get barcode 2 cluster table ---------------------------------------------
anno_tab <- seurat_object@meta.data[seurat_object@meta.data$orig.ident == sample_id,]
anno_tab$barcode <- rownames(anno_tab)
anno_tab <- anno_tab %>%
  mutate(barcode_clean = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]) %>%
  mutate(group = paste0("C", seurat_clusters)) %>%
  select(barcode_clean, group)
write.table(x = anno_tab, file = paste0(dir_out, "annotation_file.", run_id, ".txt"), quote = F, sep = "\t", row.names = F, col.names = F)
rm(seurat_object)

# input raw read count matrix ---------------------------------------------
# get matrix directory
dir_cellranger_output <- "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/Cell_Ranger/outputs/"
dir_raw <- paste0(dir_cellranger_output, sample_id, "/outs/raw_feature_bc_matrix/")

# get direct paths to data
path_matrix <- paste0(dir_raw, "matrix.mtx.gz")
path_barcodes <- paste0(dir_raw, "barcodes.tsv.gz")
path_features <- paste0(dir_raw, "features.tsv.gz")

# read in matrix
raw_exp_mat <- readMM(file = path_matrix)
feature.names = read.delim(path_features, header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(path_barcodes, header = FALSE,stringsAsFactors = FALSE)

colnames(raw_exp_mat) = barcode.names$V1
rownames(raw_exp_mat) = feature.names$V1

# expression data - remove gene not mapped in the gene position file
missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% tab_gene_order_ensembl$V1)]
missing_sym
## only missing ~750 genes, romove them
clean_exp_mat <- raw_exp_mat[!(rownames(raw_exp_mat) %in% missing_sym),]
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

