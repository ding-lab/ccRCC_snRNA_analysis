# Yige Wu @WashU Oct 2020

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
## input known TF relationship
tf2target_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/complete_tf_target_table/20201027.v1/omnipathdb.transcriptional.plus_manual.20201027.v1.tsv")
## input the DEG-TF matrix
deg2tf_wide_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/All_tumor-up_TF_DEG_withMotifs_inPromoterDACR.20201227.tsv")

# melt the DEG-TF table ---------------------------------------------------
deg2tf_long_df <- melt(deg2tf_wide_df)
deg2tf_long_filtered_df <- deg2tf_long_df %>%
  filter(value > 0) %>%
  dplyr::rename(genesymbol_deg = V1) %>%
  dplyr::rename(motifname = variable)
## format the motif names to TF names
### duplicate rows with multiple gene symbol
idx_rep <- sapply(X = deg2tf_long_filtered_df$motifname, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "(var.2)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  length_genes <- length(vec_genes)
  return(length_genes)
})
deg2tf_annotated_df <- deg2tf_long_filtered_df[rep(1:nrow(deg2tf_long_filtered_df), idx_rep),]
### get the gene symbols of the TFs
genesymbols_tf <- sapply(X = deg2tf_long_filtered_df$motifname, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "\\(var.2\\)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  return(vec_genes)
})
deg2tf_annotated_df$genesymbol_tf <- unlist(genesymbols_tf)

# merge with known relations ----------------------------------------------
deg2tf_annotated_df <- merge(x = deg2tf_annotated_df, y = tf2target_df,
                             by.x = c("genesymbol_tf", "genesymbol_deg"), by.y = c("source_genesymbol", "target_genesymbol"),
                             all.x = T)
deg2tf_annotated_df <- deg2tf_annotated_df %>%
  arrange(is_directed, desc(value), dorothea_level)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs_with_TFs_inDARs", ".Long", ".tsv")
write.table(x = deg2tf_long_filtered_df, file = file2write, quote = F, sep = "\t", row.names = )
file2write <- paste0(dir_out, "DEGs_with_TFs_inDARs", ".GeneSymbol_KnownTFTarget_Annotated", ".tsv")
write.table(x = deg2tf_annotated_df, file = file2write, quote = F, sep = "\t", row.names = )
