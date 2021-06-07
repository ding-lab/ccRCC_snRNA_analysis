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
peak2gene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/filter_peaks/filter_cnv_affected_ccRCC_specific_DACRs/20210603.v1/ccRCC_specific.DACRs.PotentialCNVEffectFiltered.20210603.v1.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/summarize_tumor_vs_pt_DEGs/20210429.v1/Tumor_DEGs.EnoughDataPoints.Consistent.20210429.v1.tsv")

# filter peaks ------------------------------------------------------------
unique(peak2gene_df$Gene) %>% length()
unique(peak2gene_df$Gene[peak2gene_df$DAP_type == "Promoter"]) %>% length()
nrow(peak2gene_df[peak2gene_df$DAP_type == "Promoter",])

degs_filter <- degs_df$genesymbol_deg[degs_df$Tumor_vs_PT == "Up"]
peaks_filtered_df <- peak2gene_df %>%
  filter(DAP_type == "Promoter") %>%
  filter(Gene %in% degs_filter) %>%
  unique()
nrow(peaks_filtered_df)
unique(peaks_filtered_df$Gene) %>% length()

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "DEG_associated_Promoter_Peaks.", run_id, ".tsv")
write.table(x = peaks_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
