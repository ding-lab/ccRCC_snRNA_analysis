# Yige Wu @WashU Apr 2021

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
# deg_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/summarize_tumor_vs_pt_DEGs/20210429.v1/Tumor_DEGs.EnoughDataPoints.Consistent.20210429.v1.tsv")
deg_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/summarize_tumor_vs_pt_DEGs/20210429.v1/Tumor_DEGs.EnoughDataPoints.20210429.v1.tsv")
deg_bulkRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/run_deg_analysis/run_tumor_vs_NAT_deg_on_cptac_ccRCC_discovery_cases/20210419.v1/Tumor_vs_NAT.glmQLFTest.OutputTables.tsv")
deg_bulkprotein_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/compare_bulk_protein_tumor_vs_normal/20200728.v1/Bulk_Protein_Tumor_vs_Normal.Wilcox.20200728.v1.tsv")

# overlap -----------------------------------------------------------------
## preprocess the bulk DEGs
deg_bulkRNA_df <- deg_bulkRNA_df %>%
  # filter(abs(logFC) > 0.5) %>%
  mutate(direction.bulkRNA = ifelse(logFC > 0, "Up", "Down")) %>%
  rename(logFC.bulkRNA = logFC) %>%
  rename(logCPM.bulkRNA = logCPM) %>%
  rename(FDR.bulkRNA = FDR) %>%
  rename(F.bulkRNA = F) %>%
  rename(PValue.bulkRNA = PValue)
deg_merged_df <- merge(x = deg_snRNA_df %>%
                          rename(direction.snRNA = Tumor_vs_PT), 
                       y = deg_bulkRNA_df, by.x = c("genesymbol_deg"), by.y = c("hgnc_symbol"), all.x = T)
## merge protein data
deg_bulkprotein_df <- deg_bulkprotein_df %>%
  mutate(direction.bulkpro = ifelse(meddiff_exp > 0, "Up", "Down")) %>%
  rename(PValue.bulkpro = p_val) %>%
  rename(meddiff_exp.bulkpro = meddiff_exp) %>%
  rename(number_values.bulkpro = number_values) %>%
  rename(FDR.bulkpro = fdr)
deg_merged_df <- merge(x = deg_merged_df, y = deg_bulkprotein_df, by.x = c("genesymbol_deg"), by.y = c("gene_symbol"), all.x = T)
deg_merged_df <- deg_merged_df %>%
  arrange(desc(Num_sig_up), desc(mean_avg_logFC))
nrow(deg_merged_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Tumor_vs_PT_DEGs.Overlap.", "snRNA.bulkRNA.Protein.", run_id, ".tsv")
write.table(x = deg_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)

