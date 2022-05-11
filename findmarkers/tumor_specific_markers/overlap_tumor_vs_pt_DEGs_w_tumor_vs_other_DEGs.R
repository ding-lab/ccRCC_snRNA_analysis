# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
degs_ccRCC_vs_pt_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_bulkRNA_protein_DEGs/20210824.v1/Tumor_vs_PT_DEGs.United.snRNA.bulkRNA.Protein.20210824.v1.tsv")
## input tumor-specific DEGs with surface annotation
# degs_ccRCC_vs_others_df
degs_ccRCC_vs_others_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/tumor_specific_markers/20210701.v1/DE_genes_filtered_surface_3DB_GTEX_HPA.txt")
## input druggable genes
druggenes_df1 <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/approved_target_ids_all.csv")
druggenes_df2 <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/01-Jul-2021-GeneSummaries.tsv")

# filter ------------------------------------------------------------------
markers_merged_df <- merge(x = degs_ccRCC_vs_others_df, 
                                           y = degs_ccRCC_vs_pt_df %>%
                                             rename(pct.allTumorcells = pct.1.allTumorcellsvsPT) %>%
                                             rename(pct.allPTcells = pct.2.allTumorcellsvsPT) %>%
                                             rename(Num_sig_up.allTumorcellsvsPT = Num_sig_up) %>%
                                             rename(Num_sig_down.allTumorcellsvsPT = Num_sig_down) %>%
                                             rename(log2FC.bulkRNA = logFC.bulkRNA) %>%
                                             rename(log2FC.bulkpro = meddiff_exp.bulkpro) %>%
                                             select(genesymbol_deg, avg_log2FC.allTumorcellsvsPT, pct.allTumorcells, pct.allPTcells,
                                                    Num_sig_up.allTumorcellsvsPT, Num_sig_down.allTumorcellsvsPT, log2FC.bulkRNA, FDR.bulkRNA, log2FC.bulkpro, FDR.bulkpro), 
                                           by.x = "Gene", by.y = "genesymbol_deg", all.x = T)
markers_merged_df$GTEX_FDR[markers_merged_df$GTEX_FDR == "FALSE"] <- NA; markers_merged_df$GTEX_FDR <- as.numeric(markers_merged_df$GTEX_FDR)
markers_merged_df$HPA_RNA_FDR[markers_merged_df$HPA_RNA_FDR == "FALSE"] <- NA; markers_merged_df$HPA_RNA_FDR <- as.numeric(markers_merged_df$HPA_RNA_FDR)
markers_merged_df$HPA_Protein_FDR[markers_merged_df$HPA_Protein_FDR == "FALSE"] <- NA; markers_merged_df$HPA_Protein_FDR <- as.numeric(markers_merged_df$HPA_Protein_FDR)

# annotate with druggability ----------------------------------------------
markers_merged_df <- merge(x = markers_merged_df, y = druggenes_df1 %>%
                                      select(`Gene Name`, `Drug IDs`) %>%
                                      rename(Drug_IDs = `Drug IDs`), by.x = c("Gene"), by.y = c("Gene Name"), all.x = T)
markers_merged_df <- merge(x = markers_merged_df, y = druggenes_df2 %>%
                                      select(name, description, gene_civic_url), by.x = c("Gene"), by.y = c("name"), all.x = T)

# filter ------------------------------------------------------------------
markers_merged_filtered_df1 <- markers_merged_df %>%
  filter((!is.na(HPA_Protein_FDR) & HPA_Protein_FDR < 0.05) | (!is.na(HPA_RNA_FDR) & HPA_RNA_FDR < 0.05) | (!is.na(GTEX_FDR) & GTEX_FDR < 0.05))

markers_merged_filtered_df2 <- markers_merged_filtered_df1 %>%
  filter(Num_sig_up.allTumorcellsvsPT >= 15 & Num_sig_down.allTumorcellsvsPT == 0)

markers_merged_filtered_df3 <- markers_merged_df %>%
  filter((!is.na(GTEX_FDR) & GTEX_FDR < 0.05)) %>%
  arrange(desc(avg_log2FC.allTumorcellsvsPT))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_vs_NonTumorcells.merge.ccRCC_vs_PT.", run_id, ".csv")
write.table(x = markers_merged_df, file = file2write, quote = F, sep = ",", row.names = F)
file2write <- paste0(dir_out, "ccRCC_vs_NonTumorcells.merge.ccRCC_vs_PT.TissueSpecific", run_id, ".tsv")
write.table(x = markers_merged_filtered_df1, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ccRCC_vs_NonTumorcells.ccRCC_highervs_PT.TissueSpecific.", run_id, ".csv")
write.table(x = markers_merged_filtered_df2, file = file2write, quote = F, sep = ",", row.names = F)
file2write <- paste0(dir_out, "ccRCC_vs_NonTumorcells.ccRCC_highervs_PT.GTEXTissueSpecific.", run_id, ".csv")
write.table(x = markers_merged_filtered_df3, file = file2write, quote = F, sep = ",", row.names = F)

