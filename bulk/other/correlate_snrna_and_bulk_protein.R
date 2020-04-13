# Yige Wu @WashU Apr 2020

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

# input denpendencies -----------------------------------------------------
## input snRNA averaged by cell type
avgexp_bycelltype_df <- fread(input = "./Resources/Analysis_Results/average_expression/adjust_averageexpression_for_bulk_correlation/20200411.v1/avgexp_bycelltype.long_data_frame.20200411.v1.tsv", data.table = F)
## load bulk protein data
protein_df <- fread("./Resources/Bulk_Processed_Data/Protein/CCRCC_PRO_tumor_PGDAC_MD_MAD_partID.txt", data.table = F)
## load id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## input differential gene result
deg_roc_df <- fread(input = "./Resources/Analysis_Results/findmarkers/filter_celltypemarkers_for_bulk_correlation/20200411.v1/celltypemarkers_filtered_for_bulk_correlation.20200411.v1.tsv", data.table = F)

# filter and map case id --------------------------------------------------
## only keep the original piece
avgexp_bycelltype_df <- avgexp_bycelltype_df %>%
  filter(aliquot %in% idmetadata_df$Aliquot.snRNA[idmetadata_df$Is_discovery_set == T & idmetadata_df$Sample_Type == "Tumor"])
## map case id
avgexp_bycelltype_df$id_case <- mapvalues(x = avgexp_bycelltype_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))

# get to do the test ------------------------------------------------------
genes2test <- unique(avgexp_bycelltype_df$gene_symbol)
length(genes2test)
genes2test <- intersect(genes2test, protein_df$Gene)
length(genes2test)
# celltype2test <- unique(avgexp_bycelltype_df$celltype)

# do correlation by gene --------------------------------------------------
## initiate vector to store the result
genes_vec <- NULL
celltypes_vec <- NULL
pvalues_vec <- NULL
rhos_vec <- NULL

for (gene_tmp in genes2test) {
  celltype2test <- deg_roc_df$cluster[deg_roc_df$gene == gene_tmp]
  celltype2test <- c(celltype2test, "All_Cells")
  for (celltype_tmp in celltype2test) {
    snrna_gene_df <- avgexp_bycelltype_df %>%
      filter(gene_symbol == gene_tmp) %>%
      filter(celltype == celltype_tmp)
    idcase2test <- snrna_gene_df$id_case
    snrna_tmp <- snrna_gene_df$avg_exp
    protein_tmp <- protein_df[protein_df$Gene == gene_tmp, idcase2test]
    protein_tmp <- as.numeric(unlist(protein_tmp))
    if (length(which(!is.na(snrna_tmp) & !is.na(protein_tmp))) > 5) {
      spearman_out <- cor.test(x = snrna_tmp, y = protein_tmp, method = "spearman")
      ## store test result
      genes_vec <- c(genes_vec, gene_tmp)
      celltypes_vec <- c(celltypes_vec, celltype_tmp)
      pvalues_vec <- c(pvalues_vec, spearman_out$p.value)
      rhos_vec <- c(rhos_vec, spearman_out$estimate)
    }
  }
}
fdrs_vec <- pvalues_vec
for (celltype_tmp in unique(spearman_result_df$celltype)) {
  fdrs_vec[celltypes_vec == celltype_tmp] <- p.adjust(p = fdrs_vec[celltypes_vec == celltype_tmp], method = "fdr")
}
spearman_result_df <- data.frame(gene_symbol = genes_vec, celltype = celltypes_vec, pvalue = pvalues_vec, rho = rhos_vec, fdr = fdrs_vec)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "spearman_result.", run_id, ".tsv")
write.table(x = spearman_result_df, file = file2write, quote = F, sep = "\t", row.names = F)
