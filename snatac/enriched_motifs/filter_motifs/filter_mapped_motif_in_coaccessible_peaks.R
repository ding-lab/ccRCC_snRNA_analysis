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
## input motif mapped
motifs_all_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Peak_Annotation/Motifs_matched.DEG_associated_Peaks.20210514.v1.tsv")
## input peak annotation
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/merge_peak_annotation/20210514.v1/Peak2Gene.20210514.v1.tsv")
## input known TF-target relations

# filter co-accessible peaks ----------------------------------------------
peaks_cap <- peaks_anno_df$Peak[peaks_anno_df$Is.CAP]
colnames(motifs_all_df)[9] <- "genesymbol_motif_nearest"
motifs_filtered_df <- motifs_all_df %>%
  dplyr::filter(Peak %in% peaks_cap)
## scenario 1: enhancer coacesssible peak mapped with enriched motif, associated with ccRCC-important gene
## scenario 2: enhancer coacesssible peak mapped with motif known to regulate the gene