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
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/annotate_peaks/annotate_BAP1_vs_NonMutant_daps_28samples/20211011.v1/BAP1_vs_NonMutant_DAP2Gene.EnhancerPromoter.20211011.v1.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/summarize_degs/unite_BAP1_vs_NonMutant_snRNA_bulkRNA_protein_DEGs/20210913.v1/BAP1_snRNA_DEGs.Consistent.CNVcorrected.20210913.v1.tsv")

# annotate and filter peaks ------------------------------------------------------------
## merge
peaks2degs_all_df <- merge(x = peaks_anno_df %>%
                             rename(avg_log2FC.snATAC = avg_log2FC) %>%
                             rename(FDR.snATAC.cnvcorrected = CNV_corr_p_adjust_bonf), 
                           y = degs_df %>%
                             select(genesymbol_deg, foldchange_type, Num_sig_up.snRNA, Num_sig_down.snRNA, avg_log2FC.snRNA, FDR.snRNA.cnvcorrected) %>%
                             rename(foldchange_type.snRNA = foldchange_type) %>%
                             mutate(DEG_direction = ifelse(avg_log2FC.snRNA > 0, "Up", "Down")),
                           by.x = c("Gene"), by.y = c("genesymbol_deg"), suffix = c(".snATAC", ".snRNA"), all = T)
peaks2degs_df <- peaks2degs_all_df %>%
  filter(!is.na(Num_sig_up.snRNA)) %>%
  filter(!is.na(peak2gene_type)) %>%
  arrange(desc(peak2gene_type))

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_vs_NonMutant_DAP2DEG.", run_id, ".tsv")
write.table(x = peaks2degs_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_vs_NonMutant_DAP_DEG_Merged.", run_id, ".tsv")
write.table(x = peaks2degs_all_df, file = file2write, quote = F, sep = "\t", row.names = F)
