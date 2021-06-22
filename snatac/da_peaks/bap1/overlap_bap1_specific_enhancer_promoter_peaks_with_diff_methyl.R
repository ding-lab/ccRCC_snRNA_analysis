# Yige Wu @WashU May 2021

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
## input daps
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/annotate_bap1_specific_daps/20210615.v1/BAP1_DAP2Gene.EnhancerPromoter.20210615.v1.tsv")
## input differential methylation probes
dm_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/methylation/compare_methylation_BAP1_tumor_vs_othertumors_katmai/20210621.v25Cores/Methyaltion_BAP1_vs_Others.Wilcox.20210621.v25Cores.tsv")
## input probe annotation
probes_anno_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/EPIC.hg38.manifest.tsv")
probe2rna_cor_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/methylation/correlate_methyl_probe_to_rna_doparallel_katmai/20210621.v25Cores/Methylation_RNA_Correlation.20210621.v25Cores.tsv")

# annotate ------------------------------------------------------------
## filter probes
dm_sig_df <- dm_df %>%
  filter(fdr < 0.05)
## annotate methylation probes
probes_anno_filtered_df <- probes_anno_df %>%
  filter(probeID %in% dm_sig_df$probeID)
dm_sig_df$probe2genes <- mapvalues(x = dm_sig_df$probeID, from = probes_anno_filtered_df$probeID, to = as.vector(probes_anno_filtered_df$gene_HGNC))
## split by gene
idx_rep <- sapply(dm_sig_df$probe2genes, function(genes) {
  genes_vec <- str_split(string = genes, pattern = ";")[[1]]
  len_genes <- length(genes_vec)
  return(len_genes)
})
genes_vec <- sapply(dm_sig_df$probe2genes, function(genes) {
  genes_vec <- str_split(string = genes, pattern = ";")[[1]]
  return(genes_vec)
})
dm2gene_df <- dm_sig_df[rep(1:nrow(dm_sig_df), idx_rep),]
dm2gene_df$probe2gene <- unlist(genes_vec)
## merge
peaks2probes_df <- merge(x = peaks_anno_df %>%
                          dplyr::rename(avg_log2FC.snATAC = avg_log2FC) %>%
                           select(Gene, peak, Count_sig.snATAC, DAP_direction, avg_log2FC.snATAC, peak2gene_type), 
                        y = dm2gene_df %>%
                          dplyr::rename(avg_log2FC.methyl = log2FC) %>%
                          dplyr::rename(fdr.methyl = fdr) %>%
                          dplyr::select(probeID, probe2gene, avg_log2FC.methyl, fdr.methyl),
                        by.x = c("Gene"), by.y = c("probe2gene"), suffix = c(".snATAC", ".methyl"))
peaks2probes_df <- unique(peaks2probes_df)
probe2rna_cor_df$fdr <- p.adjust(p = probe2rna_cor_df$p_val, method = "fdr")
peaks2probes_cor_df <- merge(x = peaks2probes_df, y = probe2rna_cor_df, by.x = c("probeID", "Gene"), by.y = c("probeID", "gene_symbol"), all.x = T)


# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_DAP2Diff_Methylation.", run_id, ".tsv")
write.table(x = peaks2probe_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_DAP2Diff_Methylation.wGeneCorrection.", run_id, ".tsv")
write.table(x = peaks2probes_cor_df, file = file2write, quote = F, sep = "\t", row.names = F)

