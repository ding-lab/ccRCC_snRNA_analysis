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
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltype13_bysample_katmai/20210701.v1/33_aliquot_merged.avgexp.SCT.data.bycelltype13_bysample.20210701.v1.tsv")
## input bulk meta data
metadata_df <- fread("./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input survival ddata
survival_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_cptac_discovery_ccRCC_survival_time/20210621.v1/CPTAC_Discovery_ccRCC_Survival_Time20210621.v1.tsv")

# reformat data -----------------------------------------------------------
exp_colnames_df <- data.frame(colname = colnames(exp_df))
exp_colnames_df <- exp_colnames_df %>%
  mutate(id_sample_cell_group = gsub(x = colname, pattern = "SCT\\.", replacement = "")) %>%
  mutate(aliquot = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,1]) %>%
  mutate(cell_group.columnname = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,2])
exp_colnames_filtered_df <- exp_colnames_df %>%
  filter(cell_group.columnname == "Tumor.cells")
exp_data_df <- exp_df[, c(exp_colnames_filtered_df$colname)]
## rename columns
exp_colnames_filtered_df$CASE_ID <- mapvalues(x = exp_colnames_filtered_df$aliquot, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Case))
colnames(exp_data_df) <- c(exp_colnames_filtered_df$CASE_ID)
#3 make colors
colors_expgroup <- RColorBrewer::brewer.pal(n = 9, name = "YlOrBr")[c(9, 3, 6)]
# namescolors_expgroup <- c("High", "Low", "Medium")
colors_expgroup <- RColorBrewer::brewer.pal(n = 9, name = "YlOrBr")[c(9, 6, 3)]

# specify gene to test ----------------------------------------------------
gene2_process <- c("CA9", "SNAP25", "EPHA6", "SHISA9", "ABLIM3", "PHKA2", "UBE2D2")

result_list <- list()
# for (gene_test in "CP") {
for (gene2_tmp in gene2_process) {
  genes_process_tmp <- c("CP", gene2_tmp)
  exp_test_wide_df <- exp_data_df[exp_df$V1 %in% genes_process_tmp,]
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
file2write <- paste0(dir_out, "Cox_2genes.3groups.", run_id, ".tsv")
write.table(x = cox_result_df, file = file2write, quote = F, sep = "\t", row.names = F)
