# Yige Wu @WashU Jun 2021

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
# peaks_up_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DA_peaks_Tumor_vs_PT_affected_byCNV_removed.tsv")
peaks_up_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/UP_Tumor_vsPT.Filtered.CNV_corrected.Annotated.20210811.tsv")
# peaks_down_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DOWN_DA_peaks_Tumor_vs_PT_affected_byCNV_removed.tsv")
peaks_down_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DOWN_Tumor_vsPT.Filtered.CNV_corrected.Annotated.20210811")
## input peak fold changes
peak2fc_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/ccRCC_specific_Filtered_peaks_byMinPct.0.1.tsv")
## input peak-to-gene annotation
peak2gene_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Peak_Annotation/28_snATACmerged_allPeaks.Annotated.20210712.tsv")

# annotate and filter peaks ------------------------------------------------------------
## annotate peaks
peaks_anno_df <- rbind(peaks_up_df %>%
                         rename(Count_sig.snATAC = Count_up) %>%
                         rename(FDR.cnvcorrected = CNV_corr_p_adjust_bonf) %>%
                         select(peak, Count_sig.snATAC, FDR.cnvcorrected) %>%
                         mutate(DAP_direction = "Up"),
                       peaks_down_df %>%
                         rename(Count_sig.snATAC = Count_down) %>%
                         rename(FDR.cnvcorrected = CNV_corr_p_adjust_bonf) %>%
                         select(peak, Count_sig.snATAC, FDR.cnvcorrected) %>%
                         mutate(DAP_direction = "Down"))
peaks_anno_df$avg_log2FC <- mapvalues(x = peaks_anno_df$peak, from = peak2fc_df$peak, to = as.vector(peak2fc_df$avg_log2FC))
peaks_anno_df$avg_log2FC[peaks_anno_df$avg_log2FC == peaks_anno_df$peak] <- NA
peaks_anno_df$avg_log2FC <- as.numeric(peaks_anno_df$avg_log2FC)
peak2gene_filtered_df <- peak2gene_df %>%
  filter(peak %in% peaks_anno_df$peak)
peaks_anno_df$annotation <- mapvalues(x = peaks_anno_df$peak, from = peak2gene_filtered_df$peak, to = as.vector(peak2gene_filtered_df$Type))
peaks_anno_df$Gene <- mapvalues(x = peaks_anno_df$peak, from = peak2gene_filtered_df$peak, to = as.vector(peak2gene_filtered_df$Gene))
peaks_anno_df <- peaks_anno_df %>%
  mutate(peak2gene_type = str_split_fixed(string = annotation, pattern = " \\(", n = 2)[,1])

peaks_anno_df %>%
  filter(!is.na(FDR.cnvcorrected)) %>%
  filter(!is.na(avg_log2FC)) %>%
  nrow()

peaks_anno_df %>%
  filter(!is.na(FDR.cnvcorrected)) %>%
  filter(!is.na(avg_log2FC)) %>%
  filter(DAP_direction = "Up") %>%
  nrow()

peaks_promoter_df <- peaks_anno_df %>%
  filter(peak2gene_type == "Promoter")

peaks_anno_df %>%
  filter(peak2gene_type == "Promoter") %>%
  nrow()

peaks_anno_df %>%
  filter(peak2gene_type != "Promoter") %>%
  nrow()

# write outupt ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_vs_PT_DAPs.Annotated.NearestGene.", run_id, ".tsv")
write.table(file = file2write, x = peaks_anno_df, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ccRCC_vs_PTDAPs.Annotated.Promoter.", run_id, ".tsv")
write.table(file = file2write, x = peak2gene_enh_pro_df, quote = F, sep = "\t", row.names = F)


