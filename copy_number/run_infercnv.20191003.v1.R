# Yige Wu @WashU Oct 2019
## for running inferCNV using integrated seruat object (raw count too big)

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# set run id ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input feature gene symbol to ensembl mapping ----------------------------
feature.names = read.delim(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/inferCNV/inputs/features.tsv.gz", header = FALSE,stringsAsFactors = FALSE)

# input gene order file ---------------------------------------------------
tab_gene_order_complete = read.delim(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/inferCNV/inputs/gencode_v21_gen_pos.complete.txt", header = FALSE,stringsAsFactors = FALSE)
tab_gene_order_complete <- tab_gene_order_complete %>%
  mutate(gene_symbol = str_split_fixed(string = V1, pattern = "\\|", n = 2)[,1]) %>%
  mutate(ensembl_id_dot = str_split_fixed(string = V1, pattern = "\\|", n = 2)[,2]) %>%
  mutate(ensembl_id = str_split_fixed(string = ensembl_id_dot, pattern = "\\.", n = 2)[,1])

# format gene order file --------------------------------------------------
tab_gene_order_complete$gene_symbol_exp <- plyr::mapvalues(tab_gene_order_complete$ensembl_id, from = feature.names$V1,to = feature.names$V2)
tab_gene_order_symbol <- tab_gene_order_complete %>%
  filter(ensembl_id != gene_symbol_exp) %>%
  filter(!duplicated(gene_symbol_exp)) %>%
  select(gene_symbol_exp, V2, V3, V4)
nrow(tab_gene_order_symbol)
write.table(x = tab_gene_order_symbol, file = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/inferCNV/inputs/gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.", run_id,".txt"), quote = F, sep = "\t", col.names = F, row.names = F)


# input integrated seurat object ------------------------------------------
int_seurat_object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/integration/integrate_seurat_objects/20190927.v1/renal_integrated.20190927.v1.RDS")

# set aliquot ID -----------------------------------------------------------
sample_id <- "CPT0001260013"

# get annotation files ------------------------------------------------------
seurat_object <- Seurat::SubsetData(int_seurat_object, idents = sample_id)

# get barcode 2 cluster table ---------------------------------------------
anno_tab <- seurat_object@meta.data[seurat_object@meta.data$orig.ident == sample_id,]
anno_tab$barcode <- rownames(anno_tab)
anno_tab <- anno_tab %>%
  select(barcode, seurat_clusters)
nrow(anno_tab)
write.table(x = anno_tab, file = paste0(dir_out, sample_id, "Barcode_Annotation.", sample_id, ".", run_id, ".txt"), quote = F, sep = "\t", row.names = F, col.names = F)

# get cells 2 process -----------------------------------------------------
sample_barcodes_int <- anno_tab$barcode
sample_barcodes_int

# get raw read count matrix -----------------------------------------------
raw_exp_mat <- seurat_object@assays$RNA@counts[, sample_barcodes_int]
dim(raw_exp_mat)

# reformat the gene order file to the gene symbols in expression matrix by the same ensembl id --------
exp_gene_symbols <- rownames(raw_exp_mat)
exp_ensembl_ids <- plyr::mapvalues(exp_gene_symbols, from = feature.names$V2,to = feature.names$V1)
## for 20k genes the gene symbols are the same
## for a lot of genes, psydogenes or lncRNA, the ensembl ids match, but the gene symbol does not match
## for the rest, the ensembl ids cannot be found in the expression data

# expression data - remove gene not mapped in the gene position file
missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% tab_gene_order_symbol$gene_symbol_exp)]
missing_sym

## only missing ~750 genes, romove them
## only keep the barcodes in the annotation files
clean_exp_mat <- raw_exp_mat[!(rownames(raw_exp_mat) %in% missing_sym), ]
rm(raw_exp_mat)
# write.table(x = clean_exp_mat, file = paste0(dir_out, "RNA_Count.", sample_id, ".", run_id, ".tsv"), quote = F, row.names = T, col.names = T, sep = "\t")
clean_exp_mat <- as.matrix(clean_exp_mat)


# check if the matrix is good for running ---------------------------------
rownames(clean_exp_mat) %in% tab_gene_order_symbol$gene_symbol_exp %>% all # Check if only gene with position left. Should be TRUE
rownames(clean_exp_mat) %>% duplicated %>% any # Check if duplicate still exist. Should be FALSE

# get ref group names -----------------------------------------------------
all_group_names <- as.vector(unique(anno_tab$seurat_clusters))
all_group_names
tumor_group_names <- c("1", "3", "4", "6", "7", "8")
ref_group_names <- all_group_names[!(all_group_names %in% tumor_group_names)]
ref_group_names

# run infercnv ------------------------------------------------------------
infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix=clean_exp_mat,
                                    annotations_file="./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/copy_number/run_infercnv/20191003.v1/CPT0086350004_annotation_file.20191003.v1.txt",
                                    delim="\t",
                                    gene_order_file= paste0(dir_out, "gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates", run_id,".txt"),
                                    ref_group_names= ref_group_names)


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=dir_out,
                             cluster_by_groups=T,
                             plot_steps=F,
                             mask_nonDE_genes = T)

