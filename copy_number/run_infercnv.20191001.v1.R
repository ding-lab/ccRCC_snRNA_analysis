# Yige Wu @WashU Oct 2019
## for running inferCNV using integrated seruat object (raw count too big)

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# input feature gene symbol to ensembl mapping ----------------------------
feature.names = read.delim(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/Cell_Ranger/outputs/CPT0075130004_notFACS/outs/raw_feature_bc_matrix/features.tsv.gz", header = FALSE,stringsAsFactors = FALSE)

# format gene order file ----------------------------------------------------
tab_gene_order_complete <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/inferCNV/inputs/gencode_v21_gen_pos.complete.txt")
tab_gene_order_ensembl <- tab_gene_order_complete %>%
  mutate(ensembl_gene_id_dot = str_split_fixed(string = V1, pattern = "\\|", n = 2)[,2]) %>%
  mutate(ensembl_gene_id = str_split_fixed(string = ensembl_gene_id_dot, pattern = "\\.", n = 2)[,1]) %>%
  select(ensembl_gene_id, V2, V3, V4)
table(tab_gene_order_ensembl$ensembl_gene_id) %>% sort() %>% tail()
## all ensembl gene ids are unique
# write.table(x = tab_gene_order_ensembl, file = paste0(dir_out, "gencode_v21_gen_pos.ensembl_gene_id.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

tab_gene_order_ensembl = read.delim(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/copy_number/run_infercnv/20191001.v1/gencode_v21_gen_pos.ensembl_gene_id.txt", header = FALSE,stringsAsFactors = FALSE)

# load seruat object ------------------------------------------------------
seurat_object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/scRNA/intergration/integrate_seurat_objects/20190927.v1/renal_integrated.20190927.v1.RDS")

# set parameters ----------------------------------------------------------
version_tmp <- 1
sample_ids <- c("CPT0075130004_notFACS", "CPT0086820004_notFACS", "CPT0075140002", "CPT0001260013", "CPT0086350004")
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
sample_id <- "CPT0086350004"

# get cells 2 process -----------------------------------------------------
sample_bc <- rownames(seurat_object@meta.data)[seurat_object@meta.data$orig.ident == sample_id]
sample_bc

# get raw read count matrix -----------------------------------------------
raw_exp_mat <- seurat_object@assays$RNA@counts[, sample_bc]
rownames(raw_exp_mat) <- plyr::mapvalues(rownames(raw_exp_mat), from = feature.names$V2,to = feature.names$V1)
rownames(raw_exp_mat) %>% head()

# expression data - remove gene not mapped in the gene position file
missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% tab_gene_order_ensembl$V1)]
missing_sym
## only missing ~750 genes
missing_row  <- match(missing_sym, rownames(raw_exp_mat))
raw_exp_mat_clean <- raw_exp_mat[-missing_row,]
rownames(raw_exp_mat_clean) %in% tab_gene_order_ensembl$V1 %>% all # Check if only gene with position left. Should be TRUE
rownames(raw_exp_mat_clean) %>% duplicated %>% any # Check if duplicate still exist. Should be FALSE
raw_exp_mat_clean <- as.matrix(raw_exp_mat_clean)
head(raw_exp_mat_clean)

# get barcode 2 cluster table ---------------------------------------------
anno_tab <- seurat_object@meta.data[seurat_object@meta.data$orig.ident == sample_id,]
anno_tab$barcode <- rownames(anno_tab)
anno_tab <- anno_tab %>%
  select(barcode, seurat_clusters)
anno_tab$seurat_clusters <- as.character(anno_tab$seurat_clusters)
write.table(x = anno_tab, file = paste0(dir_out, "annotation_file.", run_id, ".txt"), quote = F, sep = "\t", row.names = F, col.names = F)

# get ref group names -----------------------------------------------------
all_group_names <- as.vector(unique(anno_tab$seurat_clusters))
all_group_names
tumor_group_names <- c("1", "3", "4", "6", "7", "8")
ref_group_names <- all_group_names[!(all_group_names %in% tumor_group_names)]
ref_group_names

# run infercnv ------------------------------------------------------------
infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix=raw_exp_mat_clean,
                                    annotations_file="./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/copy_number/run_infercnv/20191001.v1/annotation_file.20191001.v1.txt",
                                    delim="\t",
                                    gene_order_file="./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/copy_number/run_infercnv/20191001.v1/gencode_v21_gen_pos.ensembl_gene_id.txt",
                                    ref_group_names= ref_group_names)

