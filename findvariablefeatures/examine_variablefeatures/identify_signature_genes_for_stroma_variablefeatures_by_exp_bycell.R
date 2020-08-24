# Yige Wu @WashU Aug 2020
## reference: http://jokergoo.github.io/supplementary/ComplexHeatmap-supplementary1-4/supplS2_scRNASeq/supplS2_scRNAseq.html
## Signature genes are defined as a list of genes where each gene correlates to more than 20 genes with an absolute correlation larger than 0.5.

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

# input dependencies ------------------------------------------------------
## input the variable gene list
variable_genes_df <- fread(input = "./Resources/Analysis_Results/findvariablefeatures/findvariablefeatures_cellgroup_stroma/20200811.v1/findvariablefeatures.cellgroup.Stroma.20200811.v1.tsv", data.table = F)
## input the average expression calculated (RNA)
exp_df <- fread(input = "./Resources/Analysis_Results/fetch_data/extract_exp/extract_stromacells_rnaassay_normalized_data_matrix_bycell_bystromavariablegenes/20200824.v1/stromacells_rnaassay_normalizeddata_bycell_bystromavariablegenes.20200824.v1.tsv", data.table = F)
## specify the cutoffs
cor_cutoff <- 0.5
n_cutoff <- 10

# make expression matrix --------------------------------------------------
exp_df <- exp_df %>%
  rename(gene = V1)
exp_mat <- as.matrix(exp_df[,-1])
rownames(exp_mat) <- exp_df$gene

# # Genes that are not expressed in more than half of the cells are  --------
# exp_mat = exp_mat[apply(exp_mat, 1, function(x) sum(x > 0)/length(x) > 0.5), , drop = FALSE]
# dim(exp_mat)

# identify signature genes ------------------------------------------------
## calculate gene to gene correlation
dt = cor(t(exp_mat), method = "spearman")
dt %>% head()
diag(dt) = 0
dt[abs(dt) < cor_cutoff] = 0
dt[dt < 0] = -1
dt[dt > 0] = 1
i = colSums(abs(dt)) > n_cutoff
exp_signature_mat = exp_mat[i, ,drop = FALSE]
dim(exp_signature_mat)
rownames(exp_signature_mat)
head(exp_signature_mat[1:5, 1:5])

# scale across cells -------------------------------------------------------------------
scaledexp_signature_mat = t(apply(exp_signature_mat, 1, function(x) {
  q10 = quantile(x, 0.1)
  q90 = quantile(x, 0.9)
  x[x < q10] = q10
  x[x > q90] = q90
  scale(x)
}))
head(scaledexp_signature_mat[1:5, 1:5])
colnames(scaledexp_signature_mat) = colnames(exp_signature_mat)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "stromacells_rnaassay_normalizeddata_bycell_bystromasignaturegenes.", run_id, ".tsv")
write.table(x = exp_signature_mat, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "stromacells_rnaassay_normalizeddata_bycell_bystromasignaturegenes.scaled.", run_id, ".tsv")
write.table(x = scaledexp_signature_mat, file = file2write, quote = F, sep = "\t", row.names = F)
