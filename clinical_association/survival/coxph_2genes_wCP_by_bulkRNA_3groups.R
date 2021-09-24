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
gene2_process <- c("CA9", "SNAP25", "EPHA6", "SHISA9", "ABLIM3", "PHKA2", "UBE2D2")

result_list <- list()
# for (gene_test in "CP") {
for (gene2_tmp in gene2_process) {
  genes_process_tmp <- c("CP", gene2_tmp)
  exp_test_wide_df <- exp_data_df[exp_df$gene_name %in% genes_process_tmp,]
  testdata_df <- data.frame(t(exp_test_wide_df))
  colnames(testdata_df) <- c("Expression.gene1", "Expression.gene2")
  testdata_df$CASE_ID <- colnames(exp_test_wide_df)
  testdata_df <- merge(x = testdata_df, y = survival_df, by = c("CASE_ID"), all.x = T)
  
  testdata_df <- testdata_df %>%
    filter(CASE_ID != "C3L-00359") %>%
    mutate(Expression_group.gene1 = ifelse(Expression.gene1 < quantile(x = testdata_df$Expression.gene1, probs = 0.35, na.rm = T), "Low", 
                                           ifelse(Expression.gene1 > quantile(x = testdata_df$Expression.gene1, probs = 0.65, na.rm = T), "High", "Medium"))) %>%
    mutate(Expression_group.gene2 = ifelse(Expression.gene2 < quantile(x = testdata_df$Expression.gene2, probs = 0.35, na.rm = T), "Low", 
                                           ifelse(Expression.gene2 > quantile(x = testdata_df$Expression.gene2, probs = 0.65, na.rm = T), "High", "Medium"))) %>%
    mutate(EFS_censor = (with_new_event == "Tumor Free")) %>%
    mutate(EFS = (survival_time + 9)/365) %>%
    arrange(CASE_ID, desc(Expression_group.gene1))
  testdata_df <- testdata_df[!duplicated(testdata_df$CASE_ID),]
  
  ## EFS_censor == 0 with event; == 1 without event
  ## test
  testdata_comp_df <- testdata_df %>% 
    filter(!is.na(EFS_censor) & !is.na(EFS))
  testdata_comp_df$Expression_group.gene1 <- factor(x = testdata_comp_df$Expression_group.gene1, levels = c("High", "Medium", "Low"))
  
  fit_efs <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ Expression_group.gene1 + Expression_group.gene2, data = testdata_comp_df)
  fit_efs_sum <- summary(fit_efs)
  result_list[[gene2_tmp]] <- c(genes_process_tmp,
                                fit_efs_sum$waldtest["pvalue"], fit_efs_sum$logtest["pvalue"], fit_efs_sum$sctest["pvalue"], 
                                fit_efs_sum$concordance[1],fit_efs_sum$concordance[2],
                                fit_efs_sum$coefficients[1,1], fit_efs_sum$coefficients[2,1], fit_efs_sum$coefficients[3,1], fit_efs_sum$coefficients[4,1])
}
cox_result_df <- do.call(rbind.data.frame, result_list)
colnames(cox_result_df) <- c("genesymbol1", "genesymbol2", "pvalue.wald", "pvalue.lr", "pvalue.logrank", 
                             "concordance", "se.concordance",
                             "coef.gene1.medium", "coef.gene1.low", "coef.gene2.low", "coef.gene2.medium")
cox_result_df$concordance <- as.numeric(cox_result_df$concordance)
cox_result_df$se.concordance <- as.numeric(cox_result_df$se.concordance)

cox_result_df <- cox_result_df %>%
  mutate(CI.low.concordance = concordance - 1.96*se.concordance)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cox_2genes.3groups.bulkRNA.", run_id, ".tsv")
write.table(x = cox_result_df, file = file2write, quote = F, sep = "\t", row.names = F)
