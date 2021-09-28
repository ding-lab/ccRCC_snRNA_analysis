# Yige Wu @WashU Sep 2021
## BAP1_tumorcells_vs_other_tumorcells

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
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/emt/map_motifs_to_peaks/map_motifs_to_EMT_vs_selectedEpithelial_da_peaks/20210927.v3/EMT_vs_selectedEpithelial_diff_peaks_to_motifs.20210927.v3.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_selected_EMTclusters_vs_epithelialclusters_katmai/20210924.v3/Selected_2EMTclusters_vs_5Epithelialclusters.logfc.threshold0.min.pct0.1.min.diff.pct0.AssaySCT.tsv")
## input degs
dam_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Enriched_Motifs/EMT/Score_difference.EpithelialSelectedClusters_vs_Mesenchymal.20210924.tsv")
## input known TF relationship
tf2target_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/complete_tf_target_table/20201027.v1/omnipathdb.transcriptional.plus_manual.20201027.v1.tsv")

# overlap -----------------------------------------------------------------
peaks2degs_df <- merge(x = peaks_anno_df %>%
                         filter(Type == "Promoter") %>%
                         filter(p_val_adj < 0.05) %>%
                         mutate(avg_log2FC.Epithelial_vs_Mesenchymal = avg_log2FC) %>%
                         mutate(avg_log2FC = -(avg_log2FC.Epithelial_vs_Mesenchymal)) %>%
                         select(peak, avg_log2FC, p_val_adj, Gene, Type, peak_distanceToTSS, pct.1, pct.2, motif.name, motif_coord), 
                       y = degs_df %>%
                         filter(p_val_adj < 0.05),
                       by.x = c("Gene"), by.y = c("genesymbol_deg"), suffix = c(".snATAC", ".snRNA"))
peaks2degs_df$diff_motif.mes_vs_epi <- mapvalues(x = peaks2degs_df$motif.name, from = dam_df$TF_Name, to = -as.vector(dam_df$diff))
peaks2degs_df$diff_motif.mes_vs_epi <- as.numeric(peaks2degs_df$diff_motif.mes_vs_epi)

# annotate known TF-target relations ---------------------------------------------------
## preprocess the DEG-motif table
deg2motif_df <- peaks2degs_df %>%
  mutate(TF_name = motif.name) %>%
  select(TF_name, Gene) %>%
  rename(genesymbol_deg = Gene) %>%
  unique()
### duplicate rows with multiple gene symbol
idx_rep <- sapply(X = deg2motif_df$TF_name, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "\\(var.2\\)|\\(var.3\\)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  length_genes <- length(vec_genes)
  return(length_genes)
})
deg2tf_df <- deg2motif_df[rep(1:nrow(deg2motif_df), idx_rep),]
### get the gene symbols of the TFs
genesymbols_tf <- sapply(X = deg2motif_df$TF_name, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "\\(var.2\\)|\\(var.3\\)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  return(vec_genes)
})
deg2tf_df$genesymbol_tf <- unlist(genesymbols_tf)
## preprocess the TF-targe table
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
## annotate
deg2tf_annotated_df <- merge(x = deg2tf_df, y = tf2target_uniq_df,
                             by.x = c("genesymbol_tf", "genesymbol_deg"), by.y = c("source_genesymbol", "target_genesymbol"),
                             all.x = T)

# divide by mesenchymal and epithelial ------------------------------------
## extract mesenchymal-high deg-dap-dam
dams_mes <- dam_df$TF_Name[dam_df$FDR < 0.05 & dam_df$diff < 0 & dam_df$mean_score2 > 0]
peaks2degs_mes_df <- peaks2degs_df %>%
  filter(avg_log2FC.snATAC > 0 & avg_log2FC.snRNA > 0) %>%
  filter(motif.name %in% dams_mes)
## extract epithelial-high deg-dap-dam
dams_epi <- dam_df$TF_Name[dam_df$FDR < 0.05 & dam_df$diff > 0 & dam_df$mean_score1 > 0]
peaks2degs_epi_df <- peaks2degs_df %>%
  filter(avg_log2FC.snATAC < 0 & avg_log2FC.snRNA < 0) %>%
  filter(motif.name %in% dams_epi)



