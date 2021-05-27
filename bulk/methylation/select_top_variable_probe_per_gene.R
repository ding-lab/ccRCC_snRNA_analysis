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
## input probe SD info
sd_methyl_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/methylation/calculate_gene_associated_methyl_probe_sd/20210527.v1/Methylation_Probe.SD.20210527.v1.tsv")
## input probe-to-gene info
probe2gene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/methylation/correlate_methyl_probe_to_rna/20210527.v1/Probe2Gene_HGNC.20210527.v1.tsv")
## specify the quantile cutoff
quantile_cutoff <- 0.9

# merge and select-------------------------------------------------------------------
probe2gene_sd_df <- merge(x = probe2gene_df, y = sd_methyl_df, by.x = c("probeID"), by.y = c("probe_id"), all.x = T)
probe2gene_sd_filtered_df <- probe2gene_sd_df %>%
  filter(datapoints >= 80)
sd_cutoff <- quantile(x = probe2gene_sd_filtered_df$sd, probs = quantile_cutoff)
probe2gene_sd_filtered_df2 <- probe2gene_sd_filtered_df %>%
  filter(sd >= sd_cutoff)
length(unique(probe2gene_sd_filtered_df2$probeID))
length(unique(probe2gene_sd_filtered_df2$gene_HGNC))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "top_variable_probes_per_gene.", "quantitle", quantile_cutoff, ".tsv")
write.table(x = probe2gene_sd_filtered_df2, file = file2write, quote = F, sep = "\t", row.names = F)
