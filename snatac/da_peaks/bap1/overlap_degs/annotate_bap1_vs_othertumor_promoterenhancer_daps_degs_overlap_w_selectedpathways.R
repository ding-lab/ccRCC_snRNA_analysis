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
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/annotate_bap1_specific_daps/20210615.v1/BAP1_DAP2Gene.EnhancerPromoter.20210615.v1.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/unite_BAP1_snRNA_bulkRNA_protein_DEGs/20210610.v1/BAP1_snRNA_DEGs.Consistent.CNVcorrected.20210610.v1.tsv")
deg2foldchange_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/format_all_BAP1_tumorcells_vs_other_tumorcells_CNV_corrected_degs/20210610.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.CNVcorrected.tsv")
## input selected pathways
gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pathway/unite_PBRM1_BAP1_vs_NonMutants_DAP_ORA/20210628.v2/PBRM1_BAP1_vs_NonMutants.DAPGene2TopPathway.20210628.v2.tsv")

# annotate and filter peaks ------------------------------------------------------------
## annotate DEGs
degs_df$avg_log2FC <- mapvalues(x = degs_df$genesymbol_deg, from = deg2foldchange_df$genesymbol_deg, to = as.vector(deg2foldchange_df$avg_log2FC))
degs_df$avg_log2FC <- as.numeric(degs_df$avg_log2FC)
## merge
peaks2degs_df <- merge(x = peaks_anno_df, 
                       y = degs_df %>%
                         dplyr::select(genesymbol_deg, BAP1_vs_OtherTumor_snRNA, Num_sig_up, Num_sig_down, avg_log2FC),
                       by.x = c("Gene"), by.y = c("genesymbol_deg"), suffix = c(".snATAC", ".snRNA"), all.x = T)

# add pathway -------------------------------------------------------------
pathway_dap2deg_df <- merge(x = peaks2degs_df, y = gene2pathway_df, by = c("Gene"), by.y = c("GeneSymbol"))

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_vs_othertumor_DAP2DEG.SelectedPathway.", run_id, ".tsv")
write.table(x = pathway_dap2deg_df, file = file2write, quote = F, sep = "\t", row.names = F)

## with both accessibility and expression change
## DDIT4 down in mTOR signaling, DDIT4 inhibits mTORC
## CKB, DLC1, SOS2, TRIP10, UACA down in Rho GTPase cycle
## choose DLC1 because it is reported as a TSG: https://www.nature.com/articles/s41374-018-0062-3
## COL6A2 up in focal adhession

## with only accessibility change
### TNFAIP8 up in HALLMARK_TNFA_SIGNALING_VIA_NFKB pathway
### 
peaks_anno_df %>%
  filter(!is.na(avg_log2FC)) %>%
  select(peak2gene_type, DAP_direction) %>%
  table()
# DAP_direction
# peak2gene_type Down   Up
# Enhancer 3271  344
# Promoter 2970   65

peaks2degs_df %>%
  filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA)) %>%
  select(peak2gene_type, DAP_direction) %>%
  table()
