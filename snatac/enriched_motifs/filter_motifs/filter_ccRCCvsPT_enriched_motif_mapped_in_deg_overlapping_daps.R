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
colnames(motifs_all_df)[9] <- "genesymbol_motif_nearest"
## input peak annotation
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_ccrcc_vs_pt_promoter_daps/20211004.v1/ccRCC_vs_PTDAPs.Annotated.Promoter.20211004.v1.tsv")
## specify the motifs to filter
motifs_plot <- c("NFKB2", "NFKB1", "HIF1A", "ARNT::HIF1A", "RBPJ", "MXI1", "ZNF75D", "HSF2", "NEUROD1", "SREBF2", "NEUROG2(var.2)", "TBXT", "REL", "RELA", "KLF9")

# filter co-accessible peaks ----------------------------------------------
motifs_filtered_df <- motifs_all_df %>%
  filter(Peak %in% peaks_anno_df$peak) %>%
  filter(TF_name %in% motifs_plot) %>%
  filter(Type == "Promoter")
# write outupt ------------------------------------------------------------
file2write <- paste0(dir_out, "EnrichedMotifs2PrioritizedPeaks.tsv")
write.table(x = motifs_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
