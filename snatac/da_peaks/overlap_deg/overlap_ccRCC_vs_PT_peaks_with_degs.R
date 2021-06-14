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
peaks_up_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/UP_BAP1_specific_peaks.Filtered.CNV_corrected.Annotated.20210610.tsv")
peaks_down_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/DOWN_BAP1_specific_peaks.Filtered.CNV_corrected.Annotated.20210610.tsv")
## input peak fold changes
peak2foldchange_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/BAP1_comparison_Filtered_peaks_byMinPct_MinPctDiff.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/unite_BAP1_snRNA_bulkRNA_protein_DEGs/20210610.v1/BAP1_snRNA_DEGs.Consistent.CNVcorrected.20210610.v1.tsv")
deg2foldchange_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/format_all_BAP1_tumorcells_vs_other_tumorcells_CNV_corrected_degs/20210610.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.CNVcorrected.tsv")

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
## annotate DEGs
degs_df$avg_log2FC <- mapvalues(x = degs_df$genesymbol_deg, from = deg2foldchange_df$genesymbol_deg, to = as.vector(deg2foldchange_df$avg_log2FC))
degs_df$avg_log2FC <- as.numeric(degs_df$avg_log2FC)
## merge
peaks2degs_df <- merge(x = peaks_anno_df, 
                       y = degs_df %>%
                             dplyr::select(genesymbol_deg, BAP1_vs_OtherTumor_snRNA, Num_sig_up, Num_sig_down, avg_log2FC),
                           by.x = c("Gene"), by.y = c("genesymbol_deg"), suffix = c(".snATAC", ".snRNA"))
peaks2degs_filtered_df <- peaks2degs_df %>%
  dplyr::filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA))
table(peaks2degs_filtered_df$BAP1_vs_OtherTumor_snRNA, peaks2degs_filtered_df$DAP_direction)

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_DEG_associated_Peaks.", run_id, ".tsv")
write.table(x = peaks2degs_df, file = file2write, quote = F, sep = "\t", row.names = F)
