# Yige Wu @WashU Sep 2020

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
## input tumor vs normal DA peaks
peaks2motif_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_tumor_vs_normal_da_peaks_to_celltypespecific_tfmotifs/20200916.v1/DA_peaks.Filter_for_Tumor_vs_Normal.MotifAnnotated.20200916.v1.tsv")
## input known TF relationship
tf2target_df <- fread(data.table = F, input = "./Resources/Knowledge/PPI/Transcriptional/omnipathdb.transcriptional.20200908.txt")

# annotate TF gene to DEG -------------------------------------------------
motif2deg_df <- peaks2motif_df %>%
  select(motif.name, SYMBOL, Cell_type.filename, rank_motif_by_avglogFC) %>%
  unique()
idx_rep <- sapply(X = motif2deg_df$motif.name, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "(var.2)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  length_genes <- length(vec_genes)
  return(length_genes)
})
tfgene2deg_df <- motif2deg_df[rep(1:nrow(motif2deg_df), idx_rep),]
genesymbols_tf <- sapply(X = motif2deg_df$motif.name, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "\\(var.2\\)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  return(vec_genes)
})
tfgene2deg_df$source_genesymbol <- unlist(genesymbols_tf)
motif2tfgene_df <- tfgene2deg_df %>%
  select(motif.name, source_genesymbol, Cell_type.filename, rank_motif_by_avglogFC) %>%
  unique()

tfgene2deg_df <- tfgene2deg_df %>%
  rename(target_genesymbol = SYMBOL) %>%
  select(source_genesymbol, target_genesymbol, Cell_type.filename) %>%
  unique()

# merge with known TF-target relations ------------------------------------
tfgene2deg_df <- merge(tfgene2deg_df, tf2target_df, by = c("source_genesymbol", "target_genesymbol"), all.x = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs_with_DA_peaks.", run_id, ".tsv")
write.table(x = tfgene2deg_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Motif_to_GeneSymbol.", run_id, ".tsv")
write.table(x = motif2tfgene_df, file = file2write, quote = F, sep = "\t", row.names = F)
