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
## input motif mapping
deg2motif_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/summarize_tf_motifs/count_motif_promoter_enhancer_occurance_in_up_dap_deg/20210512.v1/Motif_in_ProEnh_DAP_DEGs.20210512.v1.tsv")
calc_shared_gene_fraction <- function(vector1, vector2) {
  count_intersect <- length(which(vector1 == 1 & vector2 == 1))
  count_union <- length(which(vector1 == 1 | vector2 == 1))
  frac_shared <- count_intersect/count_union
  return(frac_shared)
}
count_shared_genes<- function(vector1, vector2) {
  count_intersect <- length(which(vector1 == 1 & vector2 == 1))
  return(count_intersect)
}
count_union_genes<- function(vector1, vector2) {
  count_union <- length(which(vector1 == 1 | vector2 == 1))
  return(count_union)
}
count_shared_genes(deg2motif_mat[, "ARNT::HIF1A"], deg2motif_mat[, "OSR1"])
count_union_genes(deg2motif_mat[, "ARNT::HIF1A"], deg2motif_mat[, "OSR1"])

# calculate shared genes between motifs -----------------------------------------
deg2motif_wide_df <- dcast(data = deg2motif_df, formula = genesymbol_deg ~ TF_name)
deg2motif_mat <- as.matrix(deg2motif_wide_df[, -1])
deg2motif_mat[deg2motif_mat > 0] <- 1

sharedgene_frac_df <- NULL
sharedgene_count_df <- NULL
for (i in colnames(deg2motif_mat)) {
  frac_i <- NULL
  count_i <- NULL
  for (j in colnames(deg2motif_mat)) {
    tmp <- calc_shared_gene_fraction(vector1 = deg2motif_mat[,i], vector2 = deg2motif_mat[,j])
    frac_i <- c(frac_i, tmp)
    tmp <- count_shared_genes(vector1 = deg2motif_mat[,i], vector2 = deg2motif_mat[,j])
    count_i <- c(count_i, tmp)
  }
  sharedgene_frac_df <- cbind(sharedgene_frac_df, frac_i)
  sharedgene_count_df <- cbind(sharedgene_count_df, count_i)
}
sharedgene_frac_df <- as.data.frame(sharedgene_frac_df); sharedgene_count_df <- as.data.frame(sharedgene_count_df)
colnames(sharedgene_frac_df) <- colnames(deg2motif_mat); rownames(sharedgene_frac_df) <- colnames(deg2motif_mat)
colnames(sharedgene_count_df) <- colnames(deg2motif_mat); rownames(sharedgene_count_df) <- colnames(deg2motif_mat)
sharedgene_frac_df <- cbind(data.frame(TF_name = rownames(sharedgene_frac_df)), sharedgene_frac_df)
sharedgene_count_df <- cbind(data.frame(TF_name = rownames(sharedgene_count_df)), sharedgene_count_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "SharedGeneFraction.tsv")
write.table(x = sharedgene_frac_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "SharedGeneCount.tsv")
write.table(x = sharedgene_count_df, file = file2write, sep = "\t", row.names = F, quote = F)
