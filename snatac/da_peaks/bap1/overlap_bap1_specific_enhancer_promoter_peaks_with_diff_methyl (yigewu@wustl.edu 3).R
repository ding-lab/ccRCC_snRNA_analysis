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
## 8512 probes that were significantly differentially methylated between BAP1 and other tumors
unique(dm2gene_df$probeID) %>% length()
unique(dm2gene_df$probe2gene) %>% length()

# process probe-to-expression correaltion ---------------------------------
probe2rna_cor_df$fdr <- p.adjust(p = probe2rna_cor_df$p_val, method = "fdr")
dm2rna_cor_df <- merge(x = dm2gene_df, y = probe2rna_cor_df, by.x = c("probeID", "probe2gene"), by.y = c("probeID", "gene_symbol"), all.x = T)
dm2rna_cor_df %>%
  select(probeID) %>%
  unique() %>%
  nrow()
dm2rna_cor_df %>%
  filter(!is.na(fdr.y)) %>%
  select(probeID) %>%
  unique() %>%
  nrow()
dm2rna_cor_df %>%
  filter(!is.na(fdr.y) & fdr.y < 0.05) %>%
  filter(rho < 0) %>%
  select(probeID) %>%
  unique() %>%
  nrow()
dm2rna_cor_df %>%
  filter(!is.na(fdr.y) & fdr.y < 0.05) %>%
  filter(rho < 0) %>%
  select(probe2gene) %>%
  unique() %>%
  nrow()

# merge peak with probe ---------------------------------------------------
peaks_anno_df %>%
  filter(peak2gene_type == "Promoter" & !is.na(avg_log2FC)) %>%
  filter(avg_log2FC < 0) %>%
  select(peak) %>%
  unique() %>%
  nrow()
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
peaks2probes_cor_df <- merge(x = peaks2probes_df, y = probe2rna_cor_df, by.x = c("probeID", "Gene"), by.y = c("probeID", "gene_symbol"), all.x = T)
peaks2probes_cor_df %>%
  filter(peak2gene_type == "Promoter" & !is.na(avg_log2FC.snATAC)) %>%
  filter(!is.na(fdr.methyl) & fdr.methyl < 0.05) %>%
  select(probeID) %>%
  unique() %>%
  nrow()
peaks2probes_cor_df %>%
  filter(peak2gene_type == "Promoter" & !is.na(avg_log2FC.snATAC)) %>%
  filter(avg_log2FC.snATAC < 0) %>%
  filter(!is.na(fdr.methyl) & fdr.methyl < 0.05) %>%
  filter(avg_log2FC.methyl > 0) %>%
  select(peak) %>%
  unique() %>%
  nrow()

peaks2probes_cor_df %>%
  filter(peak2gene_type == "Promoter" & !is.na(avg_log2FC.snATAC)) %>%
  filter(!is.na(fdr.methyl) & fdr.methyl < 0.05) %>%
  filter(!is.na(fdr.methyl) & fdr < 0.05) %>%
  filter(rho < 0) %>%
  select(probeID) %>%
  unique() %>%
  nrow()

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_DAP2Diff_Methylation.", run_id, ".tsv")
write.table(x = peaks2probe_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_DAP2Diff_Methylation.wGeneCorrection.", run_id, ".tsv")
write.table(x = peaks2probes_cor_df, file = file2write, quote = F, sep = "\t", row.names = F)

