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
dap_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/EMT/da_peaks_EpithelialSelectedClusters_vs_Mesenchymal.min.pct0.1.min.diff.pct0.logfc.threshold0.20210924.tsv")
## input peak fold changes
peak2gene_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Peak_Annotation/28_snATACmerged_allPeaks.Annotated.20210712.tsv")

# annotate and filter peaks ------------------------------------------------------------
## annotate peaks
peaks_anno_df <- dap_df %>%
  filter(p_val_adj < 0.05)
peak2gene_filtered_df <- peak2gene_df %>%
  filter(peak %in% dap_sig_df$peak)
peaks_anno_df <- merge(x = peaks_anno_df, y = peak2gene_df, by = c("peak"), all.x = T)

peaks_anno_df %>%
  filter(Type == "Promoter") %>%
  nrow()

# write outupt ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_vs_PT_DAPs.Annotated.", run_id, ".tsv")
write.table(file = file2write, x = peaks_anno_df, quote = F, sep = "\t", row.names = F)


