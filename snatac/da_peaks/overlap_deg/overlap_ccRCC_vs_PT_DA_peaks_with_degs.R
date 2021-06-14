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
peaks_up_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DA_peaks_Tumor_vs_PT_affected_byCNV_removed.tsv")
peaks_down_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DOWN_DA_peaks_Tumor_vs_PT_affected_byCNV_removed.tsv")
## input peak annotation
peak2gene_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Peak_Annotation/All_peaks_annotated_26snATAC_merged_obj.20210607.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_individual_and_CNVcorrected_DEGs/20210608.v1/Consistent.Tumor_vs_PT_DEGs.CNVcorrected.20210608.v1.tsv")

# annotate and filter peaks ------------------------------------------------------------
## annotate peaks
peaks_anno_df <- rbind(peaks_up_df %>%
                         dplyr::rename(Count_sig.snATAC = Count_up) %>%
                         dplyr::select(peak, Count_sig.snATAC, avg_log2FC) %>%
                         mutate(DAP_direction = "Up"),
                       peaks_down_df %>%
                         dplyr::rename(Count_sig.snATAC = Count_down) %>%
                         dplyr::select(peak, Count_sig.snATAC, avg_log2FC) %>%
                         mutate(DAP_direction = "Down"))
peak2gene_filtered_df <- peak2gene_df %>%
  filter(peak %in% peaks_anno_df$peak)
peaks_anno_df$annotation <- mapvalues(x = peaks_anno_df$peak, from = peak2gene_filtered_df$peak, to = as.vector(peak2gene_filtered_df$annotation))
peaks_anno_df$Gene <- mapvalues(x = peaks_anno_df$peak, from = peak2gene_filtered_df$peak, to = as.vector(peak2gene_filtered_df$SYMBOL))
peaks_anno_df <- peaks_anno_df %>%
  mutate(DAP_type = str_split_fixed(string = annotation, pattern = " \\(", n = 2)[,1])
table(peaks_anno_df$DAP_type)
## merge
peaks2degs_df <- merge(x = peaks_anno_df, 
                       y = degs_df %>%
                         rename(avg_log2FC = avg_log2FC.allTumorcellsvsPT) %>%
                         dplyr::select(genesymbol_deg, Num_sig_up, Num_sig_down, avg_log2FC),
                       by.x = c("Gene"), by.y = c("genesymbol_deg"), suffix = c(".snATAC", ".snRNA"))
peaks2degs_filtered_df <- peaks2degs_df %>%
  dplyr::filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA))

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_vs_PT_DEG_associated_DA_Peaks.", run_id, ".tsv")
write.table(x = peaks2degs_df, file = file2write, quote = F, sep = "\t", row.names = F)
