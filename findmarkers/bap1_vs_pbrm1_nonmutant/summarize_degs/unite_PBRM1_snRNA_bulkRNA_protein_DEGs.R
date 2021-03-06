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
deg_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/summarize_PBRM1_w_DoubleMutants_vs_PBRM1_NonMutant_DEGs/20210617.v1/PBRM1_DEGs.Sig.20210617.v1.tsv")
deg_bulkRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/run_deg_analysis/run_PBRM1_vs_BAP1_NonMutant_deg_on_cptac_ccRCC_discovery_cases/20210408.v1/PBRM1_Mutated_vs_BAP1_NonMutants.glmQLFTest.OutputTables.tsv")
deg_bulkprotein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/compare_bulk_protein_PBRM1_tumor_vs_othertumors/20210617.v1/Bulk_Protein_PBRM1_vs_Others.Wilcox.20210617.v1.tsv")
deg_cnvcor_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/findmarker_LR_wCNV_all_PBRM1_tumorcells_vs_PBRM1_NonMutant_cells_on_katmai/20210615.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")
deg2fc_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/findallmarker_LR_all_PBRM1_tumorcells_vs_BAP1_NonMutant_cells_on_katmai/20210615.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")

# overlap -----------------------------------------------------------------
deg_cnvcor_df$avg_log2FC.alltumorcells <- mapvalues(x = deg_cnvcor_df$row.names, from = deg2fc_df$genesymbol_deg, to = as.vector(deg2fc_df$avg_log2FC))
deg_cnvcor_df$avg_log2FC.alltumorcells <- as.numeric(deg_cnvcor_df$avg_log2FC.alltumorcells)
table(deg_snRNA_df$PBRM1_vs_OtherTumor_snRNA)
deg_merged_df <- merge(x = deg_snRNA_df, 
                       y = deg_cnvcor_df %>%
                         rename(genesymbol_deg = row.names) %>%
                         rename(FDR.cnvcorrected = FDR) %>%
                         select(genesymbol_deg, FDR.cnvcorrected, avg_log2FC.alltumorcells), by = c("genesymbol_deg"), all = T)
## preprocess the bulk DEGs
deg_bulkRNA_df2 <- deg_bulkRNA_df %>%
  # filter(abs(logFC) > 0.5) %>%
  mutate(direction.bulkRNA = ifelse(logFC > 0, "Up", "Down")) %>%
  rename(logFC.bulkRNA = logFC) %>%
  rename(logCPM.bulkRNA = logCPM) %>%
  rename(FDR.bulkRNA = FDR) %>%
  rename(F.bulkRNA = F) %>%
  rename(PValue.bulkRNA = PValue)
deg_merged_df <- merge(x = deg_merged_df, 
                       y = deg_bulkRNA_df2, by.x = c("genesymbol_deg"), by.y = c("hgnc_symbol"), all = T)
# deg_merged_df <- deg_merged_df %>%
#   arrange(desc(Num_sig_up), desc(mean_avg_log2FC))
# nrow(deg_merged_df)
# merge protein data
deg_bulkprotein_df <- deg_bulkprotein_df %>%
  mutate(direction.bulkpro = ifelse(meddiff_exp > 0, "Up", "Down")) %>%
  rename(PValue.bulkpro = p_val) %>%
  rename(meddiff_exp.bulkpro = meddiff_exp) %>%
  rename(number_bap1tumors.bulkpro = number_bap1tumors) %>%
  rename(number_othertumors.bulkpro = number_other_tumors) %>%
  rename(FDR.bulkpro = fdr)
deg_merged_df <- merge(x = deg_merged_df, y = deg_bulkprotein_df, by.x = c("genesymbol_deg"), by.y = c("gene_symbol"), all = T)

deg_filtered_df1 <- deg_merged_df %>%
  filter(!is.na(Num_up) & !is.na(FDR.cnvcorrected))
nrow(deg_filtered_df1)
table(deg_filtered_df1$PBRM1_vs_OtherTumor_snRNA)
deg_filtered_df2 <- deg_filtered_df1 %>%
  filter((Num_up == 0 & Num_sig_down >= 6.5) | (Num_down == 0 & Num_sig_up >= 6.5)) %>%
  filter(FDR.cnvcorrected < 0.05)
nrow(deg_filtered_df2)
table(deg_filtered_df2$PBRM1_vs_OtherTumor_snRNA)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_DEGs.United.", "snRNA.bulkRNA.Protein.", run_id, ".tsv")
write.table(x = deg_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "PBRM1_snRNA_DEGs.CNVcorrected.", run_id, ".tsv")
write.table(x = deg_filtered_df1, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "PBRM1_snRNA_DEGs.Consistent.CNVcorrected.", run_id, ".tsv")
write.table(x = deg_filtered_df2, file = file2write, quote = F, sep = "\t", row.names = F)

