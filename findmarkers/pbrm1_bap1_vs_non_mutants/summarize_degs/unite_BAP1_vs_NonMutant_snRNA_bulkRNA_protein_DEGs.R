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
## input the snRNA DEGs
deg_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/summarize_degs/summarize_BAP1_vs_NonMutant_DEGs/20210913.v1/BAP1_DEGs.Sig.20210913.v1.tsv")
### input the CNV correction results
deg_cnvcor_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/findmarker_LR_wCNV_all_BAP1_tumorcells_vs_NonMutant_cells_on_katmai/20210625.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")
### input overall fold change
deg_snRNA_fc_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/findmarker_LR_all_BAP1_tumorcells_vs_NonMutant_cells_on_katmai/20210624.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")
## input bulk RNA DEGs
deg_bulkRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/run_deg_analysis/run_PBRM1_BAP1_deg_on_cptac_ccRCC_discovery_cases/20210329.v1/PBRM1_BAP1_DEGs.glmQLFTest.OutputTables.tsv")
## input bulk protein DEGs
deg_bulkprotein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/protein/compare_bulk_protein_BAP1_tumor_vs_NonMutant_tumors/20210913.v1/Bulk_Protein_BAP1_vs_NonMutant.Wilcox.20210913.v1.tsv")

# processing snRNA DEGs ---------------------------------------------------
deg_cnvcor_df$avg_log2FC.snRNA <- mapvalues(x = deg_cnvcor_df$row.names, from = deg_snRNA_fc_df$genesymbol_deg, to = as.vector(deg_snRNA_fc_df$avg_log2FC))
deg_cnvcor_df$avg_log2FC.snRNA <- as.numeric(deg_cnvcor_df$avg_log2FC.snRNA)
deg_merged_df <- merge(x = deg_snRNA_df %>%
                         rename(Num_sig_up.snRNA = Num_sig_up) %>%
                         rename(Num_sig_down.snRNA = Num_sig_down) %>%
                         rename(Num_up.snRNA = Num_up) %>%
                         rename(Num_down.snRNA = Num_down) %>%
                         select(genesymbol_deg, Num_sig_up.snRNA, Num_sig_down.snRNA, Num_up.snRNA, Num_down.snRNA), 
                       y = deg_cnvcor_df %>%
                         rename(genesymbol_deg = row.names) %>%
                         rename(FDR.snRNA.cnvcorrected = FDR) %>%
                         select(genesymbol_deg, FDR.snRNA.cnvcorrected, avg_log2FC.snRNA), by = c("genesymbol_deg"), all = T)


# preprocess the bulk RNA DEGs --------------------------------------------
deg_bulkRNA_df2 <- deg_bulkRNA_df %>%
  # filter(abs(logFC) > 0.5) %>%
  filter(comparison == "BAP1_Mutated_vs_NonMutants") %>%
  mutate(direction.bulkRNA = ifelse(logFC > 0, "Up", "Down")) %>%
  rename(logFC.bulkRNA = logFC) %>%
  rename(logCPM.bulkRNA = logCPM) %>%
  rename(FDR.bulkRNA = FDR) %>%
  rename(F.bulkRNA = F) %>%
  rename(PValue.bulkRNA = PValue) %>%
  mutate(genesymbol_deg = ifelse(hgnc_symbol != "", hgnc_symbol, clone_based_ensembl_gene))
which(deg_bulkRNA_df2$hgnc_symbol == "")
deg_merged_df <- merge(x = deg_merged_df, 
                       y = deg_bulkRNA_df2, by.x = c("genesymbol_deg"), by.y = c("genesymbol_deg"), all = T)

# process bulk protein ----------------------------------------------------
deg_bulkprotein_df <- deg_bulkprotein_df %>%
  mutate(direction.bulkpro = ifelse(meddiff_exp > 0, "Up", "Down")) %>%
  dplyr::rename(PValue.bulkpro = p_val) %>%
  dplyr::rename(meddiff_exp.bulkpro = meddiff_exp) %>%
  dplyr::rename(number_group1.bulkpro = number_group1) %>%
  dplyr::rename(number_group2.bulkpro = number_group2) %>%
  dplyr::rename(FDR.bulkpro = fdr)
deg_merged_df <- merge(x = deg_merged_df, y = deg_bulkprotein_df, by.x = c("genesymbol_deg"), by.y = c("gene_symbol"), all = T)

# filter to only snRNA DEGs -----------------------------------------------
deg_filtered_df1 <- deg_merged_df %>%
  filter(!is.na(Num_sig_up.snRNA) & !is.na(FDR.snRNA.cnvcorrected)) %>%
  mutate(foldchange_type = ifelse(FDR.snRNA.cnvcorrected < 0.05, ifelse(avg_log2FC.snRNA > 0  & Num_down.snRNA == 0 & Num_sig_up.snRNA >= 5, "consistently higher in BAP1-mutants",
                                                                        ifelse(avg_log2FC.snRNA < 0  & Num_up.snRNA == 0 & Num_sig_down.snRNA >= 5, "consistently lower in BAP1-mutants",
                                                                               ifelse(Num_down.snRNA+Num_up.snRNA >= 5, "mixed fold change directions", "<50% samples w. sig. fold changes"))), "<50% samples w. sig. fold changes"))
  
nrow(deg_filtered_df1)
deg_filtered_df2 <- deg_filtered_df1 %>%
  filter(foldchange_type %in% c("consistently higher in BAP1-mutants", "consistently lower in BAP1-mutants"))
nrow(deg_filtered_df2)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_DEGs.United.", "snRNA.bulkRNA.Protein.", run_id, ".tsv")
write.table(x = deg_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_snRNA_DEGs.CNVcorrected.", run_id, ".tsv")
write.table(x = deg_filtered_df1, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_snRNA_DEGs.Consistent.CNVcorrected.", run_id, ".tsv")
write.table(x = deg_filtered_df2, file = file2write, quote = F, sep = "\t", row.names = F)

