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
## input selected probes
gene2probe_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/methylation/select_top_variable_probe_per_gene/20210527.v1/top_variable_probes_per_gene.quantitle0.9.tsv")
## input methylation
methyl_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/DNA_methylation/CPTAC_ccRCC_discovery_tumor_methylation_betavalue_probe_level_v1.0.tsv")

# average -----------------------------------------------------------------
start_time <- Sys.time()
methyl_mean_list <- parallel::mclapply(unique(gene2probe_df$gene_HGNC), function(g, g2p_df, met_df) {
  probes_tmp <- g2p_df$probeID[g2p_df$gene_HGNC == g]
  methyl_tmp_df <- met_df[met_df$Locus %in% probes_tmp,-1]
  methyl_mean_tmp_df <- colMeans(methyl_tmp_df, na.rm = T)
  return(methyl_mean_tmp_df)
}, g2p_df = gene2probe_df, met_df = methyl_df)
methyl_mean_df <- matrix(data = unlist(methyl_mean_list), byrow = T, nrow = 5)
methyl_mean_df <- as.data.frame(methyl_mean_df)
colnames(methyl_mean_df) <- names(methyl_mean_list[[1]])
end_time <- Sys.time()
end_time - start_time
## 100 loops took Time difference of 12.96927 secs for lapply
## 100 loops took Time difference of 9.062799 secs secs for mclapply
## 23490 loops: Time difference of 30.25247 mins
methyl_mean_df2 <- cbind(data.frame(gene_HGNC = unique(gene2probe_df$gene_HGNC)), methyl_mean_df)