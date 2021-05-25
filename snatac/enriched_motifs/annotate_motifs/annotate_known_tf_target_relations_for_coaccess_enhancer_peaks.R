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
deg2motif_all_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Peak_Annotation/Motifs_matched.DEG_associated_Peaks.20210514.v1.tsv")
## input peak annotation
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/merge_peak_annotation/20210514.v1/Peak2Gene.20210514.v1.tsv")
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/merge_peak_annotation/20210514.v2/Peak2Gene.20210514.v2.tsv")

peaks_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/merge_peak_annotation/20210514.v1/Peaks.20210514.v1.tsv")

# melt the DEG-TF table ---------------------------------------------------
peaks_cap <- peaks_anno_df$Peak[peaks_anno_df$Is.CAP]
colnames(deg2motif_all_df)[9] <- "genesymbol_motif_nearest"
deg2motif_df <- deg2motif_all_df %>%
  filter(Peak_Type == "Enhancer") %>%
  filter(Peak %in% peaks_cap) %>%
  select(Peak, TF_name, Gene) %>%
  rename(genesymbol_deg = Gene) %>%
  unique()
### duplicate rows with multiple gene symbol
idx_rep <- sapply(X = deg2motif_df$TF_name, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "(var.2)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  length_genes <- length(vec_genes)
  return(length_genes)
})
deg2tf_df <- deg2motif_df[rep(1:nrow(deg2motif_df), idx_rep),]
### get the gene symbols of the TFs
genesymbols_tf <- sapply(X = deg2motif_df$TF_name, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "\\(var.2\\)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  return(vec_genes)
})
deg2tf_df$genesymbol_tf <- unlist(genesymbols_tf)

# preprocess the TF-targe table -------------------------------------------
idx_rep <- sapply(X = tf2target_df$source_genesymbol, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "_")[[1]]
  length_genes <- length(vec_genes)
  return(length_genes)
})
tf2target_uniq_df <- tf2target_df[rep(x = 1:nrow(tf2target_df), as.vector(idx_rep)),]
tf2target_uniq_df <- tf2target_uniq_df %>%
  mutate(source_genesymbols = source_genesymbol)
genesymbols_tf <- sapply(X = tf2target_df$source_genesymbol, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "_")[[1]]
  return(vec_genes)
})
tf2target_uniq_df$source_genesymbol <- unlist(genesymbols_tf)

# merge with known relations ----------------------------------------------
deg2tf_annotated_df <- merge(x = deg2tf_df, y = tf2target_uniq_df,
                             by.x = c("genesymbol_tf", "genesymbol_deg"), by.y = c("source_genesymbol", "target_genesymbol"),
                             all.x = T)
deg2motif_known_df <- deg2tf_annotated_df %>%
  filter(!is.na(is_directed)) %>%
  arrange(dorothea_level) %>%
  unique()


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs_with_TFs_inDARs", ".GeneSymbol_KnownTFTarget_Annotated", ".tsv")
write.table(x = deg2tf_known_df, file = file2write, quote = F, sep = "\t", row.names = )
