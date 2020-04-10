# Yige Wu @WashU Apr 2020
## running on local
## calculate the correlation between pathway X and fraction of tumor cells with CNV Y across tumor subclusters

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input pathway score
pathway_score_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/calculate_pathway_score/unite_pathway_scores/20200407.v1/avg_pathway_scaled_avgexp.20200407.v1.tsv", data.table = )
pathway_score_df <- unique(pathway_score_df)
## input cnv fraction by chr
frac_cnv_wide_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/summarize_cnv_fraction/estimate_fraction_of_tumorcells_with_expectedcnv_perchrregion_per_manualsubcluster_using_cnvgenes/20200407.v1/fraction_of_tumorcells.expectedCNA.by_chr_region.by_manual_subcluster.20200407.v1.tsv", data.table = F)
rownames(frac_cnv_wide_df) <- frac_cnv_wide_df$manual_cluster_name

# process by chr by pathway----------------------------------------------------------
# unique(ccrcc_cna_genes_df$chr_region)
chr_region_vec <- NULL
pathway_name_vec <- NULL
spearman_rho_vec <- NULL
spearman_pvalue_vec <- NULL
for (chr_regin_tmp in colnames(frac_cnv_wide_df[,-1])) {
  ## get the subcluster names
  subclusters_tmp <- frac_cnv_wide_df$manual_cluster_name[!is.na(frac_cnv_wide_df[,chr_regin_tmp])]
  for (path_tmp in unique(pathway_score_df$pathway_name)) {
    ## get the fraction of cnvs by subcluster name
    frac_cnv_tmp <- frac_cnv_wide_df[subclusters_tmp, chr_regin_tmp]
    ## get the pathway scores by subcluster name
    pathway_score_tmp <- pathway_score_df %>%
      filter(pathway_name == path_tmp) %>%
      filter(exp_manual_subcluster_name %in% subclusters_tmp)
    rownames(pathway_score_tmp) <- pathway_score_tmp$exp_manual_subcluster_name
    pathway_score_tmp <- pathway_score_tmp[subclusters_tmp, "pathway_score"]
    ## test
    test_out <- cor.test(x = frac_cnv_tmp, y = pathway_score_tmp, method = "spearman")
    ## store results
    chr_region_vec <- c(chr_region_vec, chr_regin_tmp)
    pathway_name_vec <- c(pathway_name_vec, path_tmp)
    spearman_rho_vec <- c(spearman_rho_vec, test_out$estimate)
    spearman_pvalue_vec <- c(spearman_pvalue_vec, test_out$p.value)
  }
}


# assemble test result and adjust to fdr ----------------------------------
spearman_result <- data.frame(chr_region = chr_region_vec,
                              pathway_name = pathway_name_vec,
                              spearman_rho = spearman_rho_vec,
                              spearman_pvalue = spearman_pvalue_vec)
spearman_result$spearman_fdr <- p.adjust(p = spearman_result$spearman_pvalue, method = "fdr")

# write output -----------------------------------------------------------
file2write <- paste0(dir_out, "spearman_cor_btw_pathwayscore_and_frac_cnv.", run_id, ".tsv")
write.table(x = spearman_result, file = file2write, quote = F, sep = "\t", row.names = F)
