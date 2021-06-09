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
## input the DEG-TF matrix
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_bulkRNA_protein_DEGs/20210601.v1/Tumor_vs_PT_DEGs.United.snRNA.bulkRNA.Protein.20210601.v1.tsv")
## specify motifs
motif2gene_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Motifs_Mapped_to_Peaks/Motifs_matched.DEG_associated_Peaks.Motif_annotation.20210517.v1.tsv")
## specify genes
genes_process <- c("HK1", "HK2", "HK3", "GPI", "PFKP", "PFKL", "PFKM", "ALDOB", "ALDOA", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGK2",
                   "PGAM1", "PGAM2", "ENO1", "ENO2", "ENO3", "PKM", "PKLR", "LDHA", "LDHB", "LDHC", "LDHD",
                   "PDK1", "PDK2", "PDK3", "PDK4", "PDP1", "PDP2", "DLAT", "DLD", "PDHA1", "PDHA2", "PDHB",
                   "ACO2", "IDH2", "OGDH", "SDHB", "SDHC", "SDHD", "FH", "MDH2",
                   "ACLY", "ACACA", "FASN",
                   "PGM1", "PGM2", "UGP2", "GYS1", "PYGL", "GBE1")
## specify peaks
peaks_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DA_peaks_Tumor_vs_PT_affected_byCNV_removed.tsv")
## input motif enrichment
dam_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/summarize_tf_motifs/summarize_ccRCC_vs_PT_motifs/20210608.v1/Summarized.ccRCC_vs_PT.Motifs.20210608.v1.tsv")
motifs_degpeaks_enriched_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/filter_motifs/filter_enriched_motifs_in_deg_associated_peaks/20210601.v1/Filtered_Enriched_Motifs_in_DEG_Associated_Peaks.FDR0.001.tsv")

# filter ------------------------------------------------------------------
motif2gene_filtered_df <- motif2gene_df %>%
  filter(Gene %in% genes_process) %>%
  filter(Peak %in% peaks_df$peak) %>%
  filter(Peak_Type == "Promoter") %>%
  filter(Motif_Type == "Promoter")
motif_allpeaks_enriched_df <- dam_df %>%
  filter(foldchange_type == "consistently higher in ccRCC")
## get unique motifs
motifs_process <- unique(motif2gene_filtered_df$motif.name)
### get the gene symbols of the TFs
genesymbols_tf <- sapply(X = motifs_process, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "\\(var.2\\)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  return(vec_genes)
})
idx_rep <- sapply(X = motifs_process, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "\\(var.2\\)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  return(length(vec_genes))
})
motif2tf_df <- data.frame(motif.name = motifs_process[rep(1:length(motifs_process), idx_rep)], tf.genesymbol = unlist(genesymbols_tf))
motif2tf_de_df <- merge(x = motif2tf_df, 
                        y = deg_df %>%
                          select(genesymbol_deg, Num_sig_up, Num_sig_down, logFC.bulkRNA, FDR.bulkRNA, direction.bulkRNA, meddiff_exp.bulkpro, FDR.bulkpro, direction.bulkpro),
                        by.x = c("tf.genesymbol"), by.y = c("genesymbol_deg"), all.x = T)
## merge with the motif-gene table
motif2gene_merged_df <- merge(x = motif2gene_filtered_df %>%
                                select(Gene, Peak, Peak_Type, Motif_Type, motif.name), 
                              y = motif2tf_de_df, by= c("motif.name"), all.x = T)
motif2gene_merged_df <- motif2gene_merged_df %>%
  arrange(Gene, desc(Num_sig_up)) %>%
  unique() %>%
  mutate(Is.DEGpeaksenriched = (motif.name %in% motifs_degpeaks_enriched_df$motif.name)) %>%
  mutate(Is.Allpeaksenriched = (motif.name %in% motif_allpeaks_enriched_df$TF_Name))
  
motif2gene_prioritized1_df <- motif2gene_merged_df %>%
  filter(Is.DEGpeaksenriched | Is.Allpeaksenriched)
motif2gene_prioritized2_df <- motif2gene_prioritized1_df %>%
  filter(Num_sig_up >= 15 | (direction.bulkRNA == "Up" & !is.na(FDR.bulkRNA) & FDR.bulkRNA < 0.05)| (direction.bulkpro == "Up" & !is.na(FDR.bulkpro) & FDR.bulkpro < 0.05))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "EnrichedTFs", ".DifferentialExpression", ".tsv")
write.table(x = deg_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
