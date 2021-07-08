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
exp_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Protein/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv")
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")
## input survival ddata
survival_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_cptac_discovery_ccRCC_survival_time/20210621.v1/CPTAC_Discovery_ccRCC_Survival_Time20210621.v1.tsv")

# reformat data -----------------------------------------------------------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma")
## center data
exp_data_df <- exp_df[, metadata_filtered_df$Specimen.Label.tumor]
exp_data_df <- (exp_data_df - exp_df$ReferenceIntensity)
## rename columns
colnames(exp_data_df) <- metadata_filtered_df$CASE_ID

# specify gene to test ----------------------------------------------------
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210702.v1/ccRCC_markers.Surface.20210702.v1.tsv")
genes_process <- genes_process_df$Gene
genes_process <- genes_process[genes_process %in% exp_df$Index]

for (gene_test in genes_process) {
  
  # make combined data and test ------------------------------------------------------
  ## filter specific protein data
  exp_test_wide_df <- exp_data_df[exp_df$Index == gene_test,]
  testdata_df <- data.frame(CASE_ID = metadata_filtered_df$CASE_ID, Expression = unlist(exp_test_wide_df))
  testdata_df <- merge(x = testdata_df, y = survival_df, by = c("CASE_ID"), all.x = T)
  cutoff_exp_low <- quantile(x = testdata_df$Expression, probs = 0.35, na.rm = T); cutoff_exp_low
  cutoff_exp_high <- quantile(x = testdata_df$Expression, probs = 0.65, na.rm = T); cutoff_exp_high
  testdata_df <- testdata_df %>%
    mutate(Expression_group = ifelse(Expression < cutoff_exp_low, "Exp Low", ifelse(Expression > cutoff_exp_high, "Exp high", "Exp medium"))) %>%
    mutate(EFS_censor = (with_new_event == "Tumor Free")) %>%
    mutate(EFS = (survival_time + 9)/365)
  ## EFS_censor == 0 with event; == 1 without event
  ## test
  testdata_comp_df <- testdata_df %>%
    filter(!is.na(EFS_censor) & !is.na(EFS) & !is.na(Expression_group))
  fit_efs <- survfit(Surv(EFS, EFS_censor == 0) ~ Expression_group, data = testdata_comp_df)
  
  
  # plot --------------------------------------------------------------------
  res <- ggsurvplot(fit_efs,
                    data = testdata_comp_df,
                    conf.int = TRUE,
                    surv.median.line = "hv", pval = TRUE,
                    legend.title = "bulk protein level",
                    legend.labs = c("High", "Low", "Medium"),
                    legend = "top",
                    xlab = "Time (years)",
                    ylab = "Overall Survival",
                    # palette = c("#D95F02", "#1B9E77"),
                    ggtheme = theme_survminer(base_size = 12,
                                              base_family = "",
                                              font.main = c(12, "plain", "black"),
                                              font.submain = c(12, "plain", "black"),
                                              font.x = c(12, "plain", "black"),
                                              font.y = c(12, "plain", "black"),
                                              font.caption = c(12, "plain", "black"),
                                              font.tickslab = c(8, "plain", "black"),
                                              legend = c("top", "bottom", "left", "right", "none"),
                                              font.legend = c(8, "plain", "black")),
                    conf.int.alpha = 0.1,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "strata") # Change line type by groups
  res$table <- res$table + theme(axis.line = element_blank())
  res$plot <- res$plot + labs(title = paste0("Survival Curves by ", gene_test, " bulk protein level"))
  file2write <- paste0(dir_out, gene_test, ".png")
  png(file2write, width = 600, height = 800, res = 150)
  print(res)
  dev.off()
}
