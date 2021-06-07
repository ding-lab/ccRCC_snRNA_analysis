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
peak2gene_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Peak_Annotation/All_peaks_annotated_26snATAC_merged_obj.20210607.tsv")
peaks_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DA_peaks_Tumor_vs_PT_affected_byCNV_removed.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/summarize_tumor_vs_pt_DEGs/20210429.v1/Tumor_DEGs.EnoughDataPoints.Consistent.20210429.v1.tsv")

# annotate and filter peaks ------------------------------------------------------------
## annotate peaks
peaks_anno_df <- merge(x = peaks_df, y = peak2gene_df %>%
                         select(peak, SYMBOL, annotation), 
                       by = c("peak"), all.x = T)
peaks_anno_df <- peaks_anno_df %>%
  mutate(DAP_type = str_split_fixed(string = annotation, pattern = " \\(", n = 2)[,1])
table(peaks_anno_df$DAP_type)
## filter by DEGs
degs_filter <- degs_df$genesymbol_deg[degs_df$Tumor_vs_PT == "Up"]
peaks_filtered_df <- peaks_anno_df %>%
  filter(DAP_type == "Promoter") %>%
  rename(Gene = SYMBOL) %>%
  filter(Gene %in% degs_filter) %>%
  unique()
nrow(peaks_filtered_df)
table(peaks_filtered_df$DAP_type)
table(peaks_filtered_df$DAP_type[peaks_filtered_df$avg_log2FC > 0.6])
table(peaks_filtered_df$DAP_type[peaks_filtered_df$avg_log2FC > 0.2])

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "DEG_associated_Promoter_Peaks.", run_id, ".tsv")
write.table(x = peaks_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
