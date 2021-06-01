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
## input motifs enriched in tumor vs pt up peaks associated with DEGs
motis_uppeaks_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Enriched_Motifs/Tumor_vs_NormalPT/Enriched_motifs_inDEG_associated_Peaks.onlyPromoters_Enhancers.20210514.v1.tsv")
## input motifs enriched in tumor vs pt down peaks
motifs_downpeaks_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Enriched_Motifs/Tumor_vs_NormalPT/Enriched_motifs_in_UP_PT_vs_Tumor.min3.20210527.tsv")
## specify significance level
fdr_thres <- 0.001
# fdr_thres <- 0.05

# adjust fdr and filter ---------------------------------------------------
motis_uppeaks_df <- motis_uppeaks_df[, 1:8]
motis_uppeaks_df$fdr <- p.adjust(p = motis_uppeaks_df$pvalue, method = "fdr")
  
motifs_downpeaks_df$fdr <- p.adjust(p = motifs_downpeaks_df$pvalue, method = "fdr")

motifs_downpeaks_sig_df <- motifs_downpeaks_df %>%
  filter(fdr < fdr_thres)
motifs_downpeaks_sig <- motifs_downpeaks_sig_df$motif.name
motis_uppeaks_filtered_df <- motis_uppeaks_df %>%
  filter(fdr < fdr_thres) %>%
  filter(!(motif.name %in% motifs_downpeaks_sig)) %>%
  arrange(desc(fold.enrichment))
motis_uppeaks_filtered_df$motif.name

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Filtered_Enriched_Motifs_in_DEG_Associated_Peaks.FDR", fdr_thres, ".tsv")
write.table(x = motis_uppeaks_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
