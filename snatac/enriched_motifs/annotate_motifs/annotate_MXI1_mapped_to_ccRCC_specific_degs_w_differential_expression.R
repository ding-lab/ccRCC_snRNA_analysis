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
## specify peaks
peaks_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DA_peaks_Tumor_vs_PT_affected_byCNV_removed.tsv")
## input motif enrichment
dam_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/summarize_tf_motifs/summarize_ccRCC_vs_PT_motifs/20210608.v1/Summarized.ccRCC_vs_PT.Motifs.20210608.v1.tsv")
motifs_degpeaks_enriched_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/filter_motifs/filter_enriched_motifs_in_deg_associated_peaks/20210601.v1/Filtered_Enriched_Motifs_in_DEG_Associated_Peaks.FDR0.001.tsv")

# filter ------------------------------------------------------------------
motif2gene_filtered_df <- motif2gene_df %>%
  filter(Peak %in% peaks_df$peak) %>%
  filter(Peak_Type == "Promoter") %>%
  filter(Motif_Type == "Promoter") %>%
  filter(motif.name == "MXI1")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "MXI1", ".Mapped2DEGPromoter", ".tsv")
write.table(x = motif2gene_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)