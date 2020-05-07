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

# input dependencies ------------------------------------------------------
## input log fold change
min.pct.wilcox <- 0.1
logfc.threshold.wilcox <- 0.25
logfc_by_manualsubcluster_df <- fread(input = paste0("./Resources/Analysis_Results/findmarkers/tumor_subclusters/findallmarkers_wilcox_tumorcells_by_manualsubcluster/20200427.v1/Tumormanualsubcluster.FindAllMarkers.Wilcox.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox, ".tsv"), data.table = F)
## input pathway score
exp_df <- fread(input = "./Resources/Analysis_Results/integration/30_aliquot_integration/averageexpression/averageexpression_tumor_cells_by_manual_subcluster/20200325.v1/averageexpression_tumor_cells_by_manual_subcluster.20200325.v1.tsv", data.table = )
exp_df <- unique(exp_df)
## input cnv fraction by chr
frac_cnv_wide_df <- fread(input = "./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/estimate_fraction_of_tumorcells_with_expectedcnv_perchrregion_per_manualsubcluster_using_cnvgenes/20200407.v1/fraction_of_tumorcells.expectedCNA.by_chr_region.by_manual_subcluster.20200407.v1.tsv", data.table = F)
rownames(frac_cnv_wide_df) <- frac_cnv_wide_df$manual_cluster_name
# filter expression by selected genes -------------------------------------
## filter DEGs
logfc_by_manualsubcluster_df <- logfc_by_manualsubcluster_df %>%
  filter(p_val_adj < 0.05)
## filter expression by DEGs
exp_df <- exp_df %>%
  filter(V1 %in% logfc_by_manualsubcluster_df$gene)
## transform the data frame
exp_long_df <- melt(data = exp_df)
## get the subcluster name
exp_long_df <- exp_long_df %>%
  mutate(name_manualsubcluster = gsub(x = variable, pattern = "RNA.", replacement = "")) %>%
  select(-variable) %>%
  rename(gene_symbol = V1) %>%
  rename(exp_value = value)
# process by chr by pathway----------------------------------------------------------
# unique(ccrcc_cna_genes_df$chr_region)
chr_region_vec <- NULL
gene_symbol_vec <- NULL
spearman_rho_vec <- NULL
spearman_pvalue_vec <- NULL
for (chr_regin_tmp in colnames(frac_cnv_wide_df[,-1])) {
  ## get the subcluster names
  subclusters_tmp <- frac_cnv_wide_df$manual_cluster_name[!is.na(frac_cnv_wide_df[,chr_regin_tmp])]
  ## get the fraction of cnvs by subcluster name
  frac_cnv_tmp <- frac_cnv_wide_df[subclusters_tmp, chr_regin_tmp]
  
  for (gene_tmp in unique(exp_long_df$gene_symbol)) {
    ## get the pathway scores by subcluster name
    exp_tmp <- exp_long_df %>%
      filter(gene_symbol == gene_tmp) %>%
      filter(name_manualsubcluster %in% subclusters_tmp)
    rownames(exp_tmp) <- exp_tmp$name_manualsubcluster
    exp_tmp <- exp_tmp[subclusters_tmp, "exp_value"]
    ## test
    test_out <- cor.test(x = frac_cnv_tmp, y = exp_tmp, method = "spearman")
    ## store results
    chr_region_vec <- c(chr_region_vec, chr_regin_tmp)
    gene_symbol_vec <- c(gene_symbol_vec, gene_tmp)
    spearman_rho_vec <- c(spearman_rho_vec, test_out$estimate)
    spearman_pvalue_vec <- c(spearman_pvalue_vec, test_out$p.value)
  }
}


# assemble test result and adjust to fdr ----------------------------------
spearman_result <- data.frame(chr_region = chr_region_vec,
                              gene_symbol = gene_symbol_vec,
                              spearman_rho = spearman_rho_vec,
                              spearman_pvalue = spearman_pvalue_vec)
spearman_result$spearman_fdr <- p.adjust(p = spearman_result$spearman_pvalue, method = "fdr")

# write output -----------------------------------------------------------
file2write <- paste0(dir_out, "spearman_cor_btw_deg_avgexp_and_frac_cnv.", run_id, ".tsv")
write.table(x = spearman_result, file = file2write, quote = F, sep = "\t", row.names = F)
