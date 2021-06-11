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
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/unite_BAP1_snRNA_bulkRNA_protein_DEGs/20210610.v1/BAP1_snRNA_DEGs.Consistent.CNVcorrected.20210610.v1.tsv")

# annotate and filter peaks ------------------------------------------------------------
peaks_anno_df <- rbind(peaks_up_df %>%
                         rename(annotation = Type) %>%
                         rename(Count_sig.snATAC = Count_up) %>%
                         select(peak, annotation, Count_sig.snATAC, Gene) %>%
                         mutate(DAP_direction = "Up"),
                       peaks_down_df %>%
                         rename(annotation = Type) %>%
                         rename(Count_sig.snATAC = Count_down) %>%
                         select(peak, annotation, Count_sig.snATAC, Gene) %>%
                         mutate(DAP_direction = "Down"))
## annotate peaks
peaks_anno_df <- peaks_anno_df %>%
  mutate(DAP_type = str_split_fixed(string = annotation, pattern = " \\(", n = 2)[,1])
table(peaks_anno_df$DAP_type)
## filter by DEGs
peaks_promoter_df <- peaks_anno_df %>%
  filter(DAP_type == "Promoter")
# > table(peaks_promoter_df$DAP_direction)
# 
# Down   Up 
# 2970  134 
peaks_filtered_df <- merge(x = peaks_promoter_df, y = degs_df %>%
                             select(genesymbol_deg, BAP1_vs_OtherTumor_snRNA, Num_sig_up, Num_sig_down),
                           by.x = c("Gene", "DAP_direction"), by.y = c("genesymbol_deg", "BAP1_vs_OtherTumor_snRNA"))
# > table(peaks_filtered_df$DAP_direction)
# 
# Down   Up 
# 64    6 
peaks_filtered_df <- peaks_filtered_df %>%
  arrange(desc(Count_sig.snATAC), desc(Num_sig_down))

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_DEG_associated_Promoter_Peaks.", run_id, ".tsv")
write.table(x = peaks_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
