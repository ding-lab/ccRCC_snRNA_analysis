# Yige Wu @WashU Mar 2021

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
ora_result_df <- NULL
gsea_result_df <- NULL
for (comparison_tmp in c("BAP1_vs_Non-Mutant_unique", "PBRM1_BAP1_vs_Non-Mutant_shared", "PBRM1_vs_Non-Mutant_unique")) {
  for (deg_direction_tmp in c("up", "down")) {
    dir_clusterprofiler <- paste0("./Resources/Analysis_Results/bulk/expression/edgeR/pathway/clusterprofiler_", comparison_tmp, "_", deg_direction_tmp, "_degs/20210312.v1/")
    path_ora <- paste0(dir_clusterprofiler, "ORA_Results.tsv")
    ora_tmp_df <- fread(data.table = F, input = path_ora)
    ora_tmp_df$comparison <- comparison_tmp
    ora_tmp_df$deg_direction <- deg_direction_tmp
    ora_result_df <- rbind(ora_tmp_df, ora_result_df)
    path_gsea <- paste0(dir_clusterprofiler, "GSEA_Results.tsv")
    gsea_tmp_df <- fread(data.table = F, input = path_gsea)
    gsea_tmp_df$comparison <- comparison_tmp
    gsea_tmp_df$deg_direction <- deg_direction_tmp
    gsea_result_df <- rbind(gsea_result_df, gsea_tmp_df)
  }
}

# filter and split and combine -------------------------------------------------
## filter
gsea_result_sig_df <- gsea_result_df %>%
  filter(p.adjust < 0.05)
ora_result_sig_df <- ora_result_df %>%
  filter(p.adjust < 0.05)
## split for GSEA
degs_gsea_list <- sapply(X = gsea_result_sig_df$core_enrichment, FUN = function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = "\\/")[[1]]
  return(genes_vec)
})
gsea_row_rep <- sapply(X = gsea_result_sig_df$core_enrichment, FUN = function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = "\\/")[[1]]
  rep_len <- length(genes_vec)
  return(rep_len)
}); gsea_row_rep <- as.vector(gsea_row_rep)
gsea_gene2pathway_df <- gsea_result_sig_df[c(rep(1:nrow(gsea_result_sig_df), gsea_row_rep)),]
gsea_gene2pathway_df <- gsea_gene2pathway_df %>%
  dplyr::select(ID, Description, setSize, comparison, deg_direction)
gsea_gene2pathway_df$gene_symbol <- unlist(degs_gsea_list)
## split for ora
degs_ora_list <- sapply(X = ora_result_sig_df$geneID, FUN = function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = "\\/")[[1]]
  return(genes_vec)
})
ora_row_rep <- sapply(X = ora_result_sig_df$geneID, FUN = function(gene_string) {
  genes_vec <- str_split(string = gene_string, pattern = "\\/")[[1]]
  rep_len <- length(genes_vec)
  return(rep_len)
}); ora_row_rep <- as.vector(ora_row_rep)
ora_gene2pathway_df <- ora_result_sig_df[c(rep(1:nrow(ora_result_sig_df), ora_row_rep)),]
ora_gene2pathway_df <- ora_gene2pathway_df %>%
  dplyr::mutate(setSize = Count) %>%
  dplyr::select(ID, Description, setSize, comparison, deg_direction)
ora_gene2pathway_df$gene_symbol <- unlist(degs_ora_list)
## combine
gene2pathway_df <- rbind(ora_gene2pathway_df, gsea_gene2pathway_df)
## remove duplicates
gene2pathway_df <- gene2pathway_df %>%
  arrange(desc(deg_direction), desc(setSize))
gene2pathway_nodup_df <- gene2pathway_df[!duplicated(gene2pathway_df$gene_symbol),]

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_PBRM1_vs_Non_Mutants.DEGs2pathway.", run_id, ".tsv")
write.table(x = gene2pathway_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_PBRM1_vs_Non_Mutants.DEGs2pathway.NoDup.", run_id, ".tsv")
write.table(x = gene2pathway_nodup_df, file = file2write, quote = F, sep = "\t", row.names = F)


