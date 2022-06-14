# Yige Wu @WashU June 2022

#  set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "survival",
  "survminer"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
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
## input scores
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/other/preprocess/calculate_tumorcell_intrinsic_inflammation_score/20220612.v1/tumorcell_intrinsic_inflammation_signature_scores.20220612.v1.tsv")
## input survival ddata
survival_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_cptac_discovery_ccRCC_survival_time/20220317.v1/CPTAC_Discovery_ccRCC_Survival_Time20220317.v1.tsv")

# preprocess ------------------------------------------------------
testdata_df <- merge(x = exp_df, y = survival_df, by.x = "case", by.y = c("CASE_ID"), all.x = T)
testdata_df <- testdata_df %>%
  mutate(Expression = tumorcellintrinsic_inflamm_score)
cutoff_exp_low <- quantile(x = testdata_df$Expression, probs = 0.25, na.rm = T); cutoff_exp_low
cutoff_exp_high <- quantile(x = testdata_df$Expression, probs = 0.75, na.rm = T); cutoff_exp_high
testdata_df <- testdata_df %>%
  mutate(Expression_group = ifelse(Expression <= cutoff_exp_low, "Low", ifelse(Expression >= cutoff_exp_high, "High", "Medium")))
table(testdata_df$Expression_group)

# test overall survival ---------------------------------------------------
## EFS_censor == 0 with event; == 1 without event
## test
testdata_comp_df <- testdata_df %>%
  mutate(surv_time = (OS_time + 9))  %>%
  mutate(surv_status = ifelse(OS_status == "censored", 1, 2)) %>%
  filter(!is.na(surv_status) & !is.na(surv_time) & !is.na(Expression_group)) %>%
  filter(Expression_group != "Medium")
fit_efs <- survfit(Surv(surv_time, surv_status) ~ Expression_group, data = testdata_comp_df)

res <- ggsurvplot(fit_efs,
                  data = testdata_comp_df,
                  conf.int = TRUE,
                  surv.median.line = "hv", pval = TRUE,
                  legend.title = paste0("tumor-cell-intrinsic\ninflammation score\n(mRNA)"),
                  legend.labs = c("High", "Low"),
                  legend = "top",
                  xlab = "Time (days)",
                  ylab = "Overall Survival",
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
file2write <- paste0(dir_out, "OS.pdf")
pdf(file2write, width = 4, height = 5, useDingbats = F)
print(res)
dev.off()

# test progression-free survival ---------------------------------------------------
## test
testdata_comp_df <- testdata_df %>%
  mutate(surv_time = (PFS_time + 9))  %>%
  mutate(surv_status = ifelse(PFS_status == "censored", 1, 2)) %>%
  filter(!is.na(surv_status) & !is.na(surv_time) & !is.na(Expression_group)) %>%
  filter(Expression_group != "Medium")
fit_efs <- survfit(Surv(surv_time, surv_status) ~ Expression_group, data = testdata_comp_df)

res <- ggsurvplot(fit_efs,
                  data = testdata_comp_df,
                  conf.int = TRUE,
                  surv.median.line = "hv", pval = TRUE,
                  legend.title = paste0("tumor-cell-intrinsic\ninflammation score\n(mRNA)"),
                  legend.labs = c("High", "Low"),
                  legend = "top",
                  xlab = "Time (days)",
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
file2write <- paste0(dir_out, "PFS.pdf")
pdf(file2write, width = 4, height = 5, useDingbats = F)
print(res)
dev.off()

