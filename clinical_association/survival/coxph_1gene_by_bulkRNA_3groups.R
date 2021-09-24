# Yige Wu @WashU Jun 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(survival)
library(survminer)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input protein data
# exp_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/RNA_rpkm_tumor_normal.tsv")
exp_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/CPTAC_ccRCC_discovery_tumor_mRNA_FPKM_UQ_log2_v1.0.tsv")
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")
## input survival ddata
survival_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_cptac_discovery_ccRCC_survival_time/20210621.v1/CPTAC_Discovery_ccRCC_Survival_Time20210621.v1.tsv")

# make combined data and test ------------------------------------------------------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma") %>%
  mutate(Sample_ID = paste0(CASE_ID, "-T"))
## subset data
exp_data_df <- exp_df[, metadata_filtered_df$Sample_ID]
## rename columns
colnames(exp_data_df) <- metadata_filtered_df$CASE_ID
#3 make colors
colors_expgroup <- RColorBrewer::brewer.pal(n = 9, name = "YlOrBr")[c(9, 3, 6)]

# specify gene to test ----------------------------------------------------
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210824.v1/ccRCC_markers.Surface.20210824.v1.tsv")
genes_process <- genes_process_df$Gene
genes_process <- genes_process[!(genes_process %in% c("DPP6", "CPNE8", "EFNA5", "MGLL", "SPIRE1", "SPIRE1", "PLCB1", "OSMR", "SORBS1", "ANO6", "EPB41", "PAM", "RHEX"))]

result_list <- list()
# for (gene_test in "CP") {
for (gene_test in genes_process) {
  ## filter specific protein data
  exp_test_wide_df <- exp_data_df[exp_df$gene_name == gene_test,]
  testdata_df <- data.frame(CASE_ID = metadata_filtered_df$CASE_ID, Expression = unlist(exp_test_wide_df))
  testdata_df <- merge(x = testdata_df, y = survival_df, by = c("CASE_ID"), all.x = T)
  cutoff_exp_low <- quantile(x = testdata_df$Expression, probs = 0.35, na.rm = T); cutoff_exp_low
  cutoff_exp_high <- quantile(x = testdata_df$Expression, probs = 0.65, na.rm = T); cutoff_exp_high
  testdata_df <- testdata_df %>%
    mutate(Expression_group = ifelse(Expression < cutoff_exp_low, "Low", ifelse(Expression > cutoff_exp_high, "High", "Medium"))) %>%
    mutate(EFS_censor = (with_new_event == "Tumor Free")) %>%
    mutate(EFS = (survival_time + 9)/365)
  ## EFS_censor == 0 with event; == 1 without event
  ## test
  testdata_comp_df <- testdata_df %>%
    filter(!is.na(EFS_censor) & !is.na(EFS) & !is.na(Expression_group))
  fit_efs <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ Expression_group, data = testdata_comp_df)
  fit_efs_sum <- summary(fit_efs)
  result_list[[gene_test]] <- c(fit_efs_sum$waldtest["pvalue"], fit_efs_sum$logtest["pvalue"], fit_efs_sum$sctest["pvalue"],
                                fit_efs_sum$coefficients[1,1], fit_efs_sum$coefficients[2,1], fit_efs_sum$concordance[1])
  
}
cox_result_df <- do.call(rbind.data.frame, result_list)
cox_result_df <- cbind(genes_process, cox_result_df)
colnames(cox_result_df) <- c("genesymbol", "pvalue.wald", "pvalue.lr", "pvalue.logrank", "coef.medium", "coef.low", "concordance")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cox_1gene.3groups.", run_id, ".tsv")
write.table(x = cox_result_df, file = file2write, quote = F, sep = "\t", row.names = F)
