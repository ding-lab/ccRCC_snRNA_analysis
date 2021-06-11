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
## input all DAPs
peaks_up_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/UP_BAP1_specific_peaks.Filtered.CNV_corrected.Annotated.20210610.tsv")
peaks_down_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/DOWN_BAP1_specific_peaks.Filtered.CNV_corrected.Annotated.20210610.tsv")
## input DAPs for the cell lines
daps_cellline_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/DA_Peaks.BAP1_vs_Control.786O_CellLines.min.diff.pct.0min.pct.0.1logfc.threshold.0.Annotated.tsv")

# annotate and filter peaks ------------------------------------------------------------
## annotate peaks
peaks_anno_df <- rbind(peaks_up_df %>%
                         dplyr::rename(annotation = Type) %>%
                         dplyr::rename(Count_sig.snATAC = Count_up) %>%
                         dplyr::select(peak, annotation, Count_sig.snATAC, Gene) %>%
                         mutate(DAP_direction = "Up"),
                       peaks_down_df %>%
                         dplyr::rename(annotation = Type) %>%
                         dplyr::rename(Count_sig.snATAC = Count_down) %>%
                         dplyr::select(peak, annotation, Count_sig.snATAC, Gene) %>%
                         mutate(DAP_direction = "Down"))
peaks_anno_df$avg_log2FC <- mapvalues(x = peaks_anno_df$peak, from = peak2foldchange_df$peak, to = as.vector(peak2foldchange_df$avg_log2FC))
peaks_anno_df$avg_log2FC <- as.numeric(as.vector(peaks_anno_df$avg_log2FC))
peaks_anno_df <- peaks_anno_df %>%
  mutate(DAP_type = str_split_fixed(string = annotation, pattern = " \\(", n = 2)[,1])
table(peaks_anno_df$DAP_type)
## merge
peaks_merged_df <- merge(x = peaks_anno_df %>%
                           dplyr::select(-annotation) %>%
                           dplyr::filter(!is.na(avg_log2FC)), 
                       y = daps_cellline_df %>%
                         dplyr::filter(p_val_adj < 0.05) %>%
                         dplyr::mutate(DAP_type = str_split_fixed(string = Type, pattern = " \\(", n = 2)[,1]) %>%
                         dplyr::rename(p_val_adj.snATAC.cellline = p_val_adj) %>%
                         dplyr::select(Gene, peak, avg_log2FC, p_val_adj.snATAC.cellline, DAP_type),
                       by = c("Gene"), suffix = c(".snATAC.tumortissue", ".snATAC.cellline"))

peaks_consistent_df <- peaks_merged_df %>%
  dplyr::filter((avg_log2FC.snATAC.tumortissue < 0 & avg_log2FC.snATAC.cellline < 0) | (avg_log2FC.snATAC.tumortissue > 0 & avg_log2FC.snATAC.cellline > 0) )
peaks_consistent_promoter_df <- peaks_consistent_df %>%
  dplyr::filter(DAP_type.snATAC.tumortissue == "Promoter" & DAP_type.snATAC.cellline == "Promoter")

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_associated_Peaks.HumanTumorTissue_CellLine.", run_id, ".tsv")
write.table(x = peaks_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Consistent.BAP1_associated_Peaks.HumanTumorTissue_CellLine.", run_id, ".tsv")
write.table(x = peaks_consistent_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Consistent.Promoter.BAP1_associated_Peaks.HumanTumorTissue_CellLine.", run_id, ".tsv")
write.table(x = peaks_consistent_promoter_df, file = file2write, quote = F, sep = "\t", row.names = F)

