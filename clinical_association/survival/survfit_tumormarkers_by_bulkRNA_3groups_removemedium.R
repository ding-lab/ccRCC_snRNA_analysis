# Yige Wu @WashU Jun 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "survival",
  "survminer"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input protein data
# exp_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/RNA_rpkm_tumor_normal.tsv")
exp_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/CPTAC_ccRCC/CPTAC_ccRCC_discovery_tumor_mRNA_FPKM_UQ_log2_v1.0.tsv")
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")
## input survival ddata
# survival_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_cptac_discovery_ccRCC_survival_time/20210621.v1/CPTAC_Discovery_ccRCC_Survival_Time20210621.v1.tsv")
survival_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_cptac_discovery_ccRCC_survival_time/20220317.v1/CPTAC_Discovery_ccRCC_Survival_Time20220317.v1.tsv")

# make combined data and test ------------------------------------------------------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma") %>%
  mutate(Sample_ID = paste0(CASE_ID, "-T"))
## subset data
exp_data_df <- exp_df[, metadata_filtered_df$Sample_ID]
## rename columns
colnames(exp_data_df) <- metadata_filtered_df$CASE_ID
#3 make colors
colors_expgroup <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[c(9, 4, 6)]

# specify gene to test ----------------------------------------------------
genes_process <- c("ABCC3", "ABLIM3", "COL23A1", "CP", "EGFR", "ENPP3", "EPHA6", "FTO", "KCTD3", "NDRG1", "PCSK6", "PHKA2", "PLEKHA1", "PLIN2", "SEMA6A", "SHISA9", "SLC6A3", "SNAP25", "TGFA", "UBE2D2" )
os_pvalue_vec <- NULL
pfs_pvalue_vec <- NULL
fontsize <- 14
for (gene_test in genes_process) {
  ## filter specific protein data
  exp_test_wide_df <- exp_data_df[exp_df$gene_name == gene_test,]
  testdata_df <- data.frame(CASE_ID = metadata_filtered_df$CASE_ID, Expression = unlist(exp_test_wide_df))
  testdata_df <- merge(x = testdata_df, y = survival_df, by = c("CASE_ID"), all.x = T)
  cutoff_exp_low <- quantile(x = testdata_df$Expression, probs = 0.35, na.rm = T); cutoff_exp_low
  cutoff_exp_high <- quantile(x = testdata_df$Expression, probs = 0.65, na.rm = T); cutoff_exp_high
  testdata_df <- testdata_df %>%
    mutate(Expression_group = ifelse(Expression < cutoff_exp_low, "Low", ifelse(Expression > cutoff_exp_high, "High", "Medium")))
  
  
  # test overall survival ---------------------------------------------------
  testdata_comp_df <- testdata_df %>%
    mutate(surv_time = (OS_time + 9)/365)  %>%
    mutate(surv_status = ifelse(OS_status == "censored", 1, 2)) %>%
    filter(!is.na(surv_status) & !is.na(surv_time) & !is.na(Expression_group)) %>%
    filter(Expression_group != "Medium")
  fit_efs <- survfit(Surv(surv_time, surv_status) ~ Expression_group, data = testdata_comp_df)
  os_pvalue_vec <- c(os_pvalue_vec, surv_pvalue(fit_efs)[1,"pval"])
  
  res <- ggsurvplot(fit_efs,
                    data = testdata_comp_df,
                    conf.int = TRUE,
                    surv.median.line = "hv", pval = TRUE,
                    legend.title = paste0(gene_test, " expression\n(bulk RNA-seq)"),
                    legend.labs = c("High", "Low"),
                    legend = "top",
                    xlab = "Time (years)",
                    ylab = "Overall Survival",
                    palette = colors_expgroup[-3],
                    ggtheme = theme_survminer(base_size = fontsize,
                                              base_family = "",
                                              font.main = c(fontsize, "plain", "black"),
                                              font.submain = c(fontsize, "plain", "black"),
                                              font.x = c(fontsize, "plain", "black"),
                                              font.y = c(fontsize, "plain", "black"),
                                              font.caption = c(fontsize, "plain", "black"),
                                              font.tickslab = c(fontsize, "plain", "black"),
                                              legend = c("top", "bottom", "left", "right", "none"),
                                              font.legend = c(fontsize, "plain", "black")),
                    conf.int.alpha = 0.1,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    risk.table.height = 0.3,
                    linetype = "strata") # Change line type by groups
  res$table <- res$table + theme(axis.line = element_blank())
  file2write <- paste0(dir_out, gene_test, ".OS.pdf")
  pdf(file2write, width = 4, height = 4.75, useDingbats = F)
  print(res)
  dev.off()
  
  # test progression-free survival ---------------------------------------------------
  ## test
  testdata_comp_df <- testdata_df %>%
    mutate(surv_time = (PFS_time + 9)/365)  %>%
    mutate(surv_status = ifelse(PFS_status == "censored", 1, 2)) %>%
    filter(!is.na(surv_status) & !is.na(surv_time) & !is.na(Expression_group)) %>%
    filter(Expression_group != "Medium")
  fit_efs <- survfit(Surv(surv_time, surv_status) ~ Expression_group, data = testdata_comp_df)
  pfs_pvalue_vec <- c(pfs_pvalue_vec, surv_pvalue(fit_efs)[1,"pval"])
  
  res <- ggsurvplot(fit_efs,
                    data = testdata_comp_df,
                    conf.int = TRUE,
                    surv.median.line = "hv", pval = TRUE,
                    legend.title = paste0(gene_test, " expression\n(snRNA-seq)"),
                    legend.labs = c("High", "Low"),
                    legend = "top",
                    xlab = "Time (years)",
                    ylab = "Progression-free Survival",
                    palette = c("#800026", "#FEB24C"),
                    ggtheme = theme_survminer(base_size = 12,
                                              base_family = "",
                                              font.main = c(12, "plain", "black"),
                                              font.submain = c(12, "plain", "black"),
                                              font.x = c(12, "plain", "black"),
                                              font.y = c(12, "plain", "black"),
                                              font.caption = c(12, "plain", "black"),
                                              font.tickslab = c(12, "plain", "black"),
                                              legend = c("top", "bottom", "left", "right", "none"),
                                              font.legend = c(12, "plain", "black")),
                    conf.int.alpha = 0.1,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "strata") # Change line type by groups
  res$table <- res$table + theme(axis.line = element_blank())
  file2write <- paste0(dir_out, gene_test, ".PFS.pdf")
  pdf(file2write, width = 4, height = 5, useDingbats = F)
  print(res)
  dev.off()
}
survfit_result_df <- data.frame(gene = genes_process, 
                                OS_pvalue = os_pvalue_vec, 
                                OS_fdr = p.adjust(p = os_pvalue_vec, method = "fdr"),
                                PFS_pvalue = pfs_pvalue_vec, 
                                PFS_fdr = p.adjust(p = pfs_pvalue_vec, method = "fdr"))


file2write <- paste0(dir_out, "surfit_results.", run_id, ".tsv")
write.table(x = survfit_result_df, file = file2write, quote = F, sep = "\t", row.names = F)

