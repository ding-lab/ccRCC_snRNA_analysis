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
## input the correlation result
cor_result_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/methylation/correlate_methyl_probe_to_rna_parallel_test/20210527.v1/Methylation_RNA_Correlation.20210527.v1.tsv")

# filter ------------------------------------------------------------------
cor_result_filtered_df <- cor_result_df %>%
  filter(p_val < 0.05) %>%
  filter(rho < 0) %>%
  arrange(gene_symbol)
unique(cor_result_filtered_df$gene_symbol)
topprobe_df <- cor_result_filtered_df %>%
  group_by(gene_symbol) %>%
  slice_min(n = 1, order_by = rho)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "TopRNACorrelated_MethylationProbes.Selected32",  ".tsv")
write.table(x = topprobe_df, file = file2write, quote = F, sep = "\t", row.names = F)
