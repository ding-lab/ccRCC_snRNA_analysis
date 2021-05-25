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
## input all dap + caps
peak2gene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/merge_peak_annotation/20210517.v1/Peak2Gene.20210517.v1.tsv")
## input peak up count
peak_count_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Differential_Peaks/20210510/da_up_peaks_Tumor_vs_PT.annotated.20210510.tsv")
## input selected pathway degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/plotting/heatmap_enriched_motif2dap_cap_deg/20210518.v1/DEGs.SelectedPathways.tsv")

# filter peaks ------------------------------------------------------------
degs_filter <- degs_df$Gene
peaks_filtered_df <- peak2gene_df %>%
  filter(Gene %in% degs_filter) %>%
  unique()
unique(peaks_filtered_df$Peak) %>% length()
peaks_filtered_df$Count_up <- mapvalues(x = peaks_filtered_df$Peak, from = peak_count_df$peak, to = as.vector(peak_count_df$Count_up))
peaks_filtered_df$Count_up[peaks_filtered_df$Count_up == peaks_filtered_df$Peak] <- 0
peaks_filtered_df$Count_up <- as.numeric(peaks_filtered_df$Count_up)
peaks_filtered_df <- peaks_filtered_df %>%
  arrange(desc(Peak_Type), desc(Count_up))

peaks_prioritized_df <- peaks_filtered_df %>%
  filter(Peak_Type == "Promoter" & Count_up > 10)

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "DEG_associated_Peaks.", run_id, ".tsv")
write.table(x = peaks_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "DEG_associated_Peaks.Promoter.Prioritized.", run_id, ".tsv")
write.table(x = peaks_prioritized_df, file = file2write, quote = F, sep = "\t", row.names = F)
