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
## input daps
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pbrm1/annotate_pbrm1_specific_daps/20210616.v1/PBRM1_DAP2Gene.EnhancerPromoter.20210616.v1.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/unite_PBRM1_snRNA_bulkRNA_protein_DEGs/20210617.v1/PBRM1_snRNA_DEGs.Consistent.CNVcorrected.20210617.v1.tsv")
deg2foldchange_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/findallmarker_LR_all_PBRM1_tumorcells_vs_BAP1_NonMutant_cells_on_katmai/20210615.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")

# annotate and filter peaks ------------------------------------------------------------
## annotate DEGs
degs_df$avg_log2FC <- mapvalues(x = degs_df$genesymbol_deg, from = deg2foldchange_df$genesymbol_deg, to = as.vector(deg2foldchange_df$avg_log2FC))
degs_df$avg_log2FC <- as.numeric(degs_df$avg_log2FC)
## merge
peaks2degs_df <- merge(x = peaks_anno_df, 
                       y = degs_df %>%
                             dplyr::select(genesymbol_deg, PBRM1_vs_OtherTumor_snRNA, Num_sig_up, Num_sig_down, avg_log2FC),
                           by.x = c("Gene"), by.y = c("genesymbol_deg"), suffix = c(".snATAC", ".snRNA"))
peaks2degs_filtered_df <- peaks2degs_df %>%
  dplyr::filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA))
table(peaks2degs_filtered_df$PBRM1_vs_OtherTumor_snRNA, peaks2degs_filtered_df$DAP_direction)
nrow(peaks2degs_filtered_df)

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_DAP2DEG.", run_id, ".tsv")
write.table(x = peaks2degs_df, file = file2write, quote = F, sep = "\t", row.names = F)
