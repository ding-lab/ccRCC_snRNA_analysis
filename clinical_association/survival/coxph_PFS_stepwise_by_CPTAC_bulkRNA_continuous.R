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
library(My.stepwise)
## set run id
version_tmp <- 5
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
survival_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_cptac_discovery_ccRCC_survival_time/20210920.v1/CPTAC_Discovery_ccRCC_Survival_Time20210920.v1.tsv")
## input marker gene
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210824.v1/ccRCC_markers.Surface.20210824.v1.tsv")

# make combined data and test ------------------------------------------------------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma") %>%
  mutate(Sample_ID = paste0(CASE_ID, "-T"))
## subset data
exp_data_df <- exp_df[, metadata_filtered_df$Sample_ID]
## rename columns
colnames(exp_data_df) <- metadata_filtered_df$CASE_ID

# specify gene to test ----------------------------------------------------
genes_process <- genes_process_df$Gene
genes_process <- genes_process[!(genes_process %in% c("DPP6", "CPNE8", "EFNA5", "MGLL", "SPIRE1", "SPIRE1", "PLCB1", "OSMR", "SORBS1", "ANO6", "EPB41", "PAM", "RHEX"))]

# pre-process -----------------------------------------------------------------
exp_test_wide_df <- exp_data_df[exp_df$gene_name %in% genes_process,]
testdata_df <- data.frame(t(exp_test_wide_df))
colnames(testdata_df) <- exp_df$gene_name[exp_df$gene_name %in% genes_process]
testdata_df$CASE_ID <- colnames(exp_test_wide_df)
testdata_df <- merge(x = testdata_df, y = survival_df, by = c("CASE_ID"), all.x = T)
testdata_df <- testdata_df %>%
  dplyr::mutate(EFS_censor = (with_new_event == "With Tumor")) %>%
  dplyr::mutate(EFS = (survival_time + 9)/365) %>%
  arrange(CASE_ID) %>% 
  filter(!is.na(EFS_censor) & !is.na(EFS))
genes_test <- exp_df$gene_name[exp_df$gene_name %in% genes_process]

# process expression only -------------------------------------------------
file2write <- paste0(dir_out, "Expressiononly.Stepwise.Cox.", run_id, ".txt")
sink(file2write)
# cat(paste0("sle: 0.15; sls: 0.15\n"))
cat(paste0("sle: 0.25; sls: 0.25\n"))
My.stepwise.coxph(Time = "EFS", Status = "EFS_censor", variable.list = genes_test, data = testdata_df, sle = 0.25, sls = 0.25)
sink()

# process expression + basic patient info -------------------------------------------------
file2write <- paste0(dir_out, "Expression_plusAgeSexetc.Stepwise.Cox.", run_id, ".txt")
sink(file2write)
My.stepwise.coxph(Time = "EFS", Status = "EFS_censor", variable.list = c("age", "sex.ismale", "stage.numeric", "grade.numeric", genes_test),  data = testdata_df, sle = 0.3, sls = 0.3)
sink()

# evaluate V1 -------------------------------------------------------------
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ CP + SHISA9, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "CP_SHISA9.Cox.", run_id, ".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ stage.numeric + grade.numeric + sex.ismale + age + CP + SHISA9, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "Stage_Grade_Sex_Age.CP_SHISA9.Cox.", run_id, ".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# # manually test best model for V2 ------------------------------------------------
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ SHISA9 + ABCC3 + NDRG1 + SEMA6A + EPHA6 + KCTD3 + UBE2D2, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "1.ExpressionOnly.DropPLIN2",".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ SHISA9 + NDRG1 + SEMA6A + EPHA6 + KCTD3 + UBE2D2, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "2.ExpressionOnly.DropPLIN2_ABCC3",".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ SHISA9 + NDRG1 + SEMA6A + EPHA6 + UBE2D2, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "3.ExpressionOnly.DropPLIN2_ABCC3_KTCD3",".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ SHISA9 + NDRG1 + SEMA6A + UBE2D2, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "4.ExpressionOnly.DropPLIN2_ABCC3_KTCD3_EPHA6",".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ SHISA9 + NDRG1 + SEMA6A, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "5.ExpressionOnly.DropPLIN2_ABCC3_KTCD3_EPHA6_UBE2D2",".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ SHISA9 + NDRG1, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "6.ExpressionOnly.DropPLIN2_ABCC3_KTCD3_EPHA6_UBE2D2_SEMA6A",".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ stage.numeric + grade.numeric + sex.ismale + age + SHISA9 + NDRG1 + stage.numeric , data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "7.ExpressionOnly.DropPLIN2_ABCC3_KTCD3_EPHA6_UBE2D2_SEMA6A",".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# ## drop PHKA2
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ stage.numeric + EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + sex.ismale + EGFR + CP + PLEKHA1, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "Drop_PHKA2.Cox.", run_id, ".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# ## drop PHKA2
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ stage.numeric + EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + sex.ismale + EGFR + CP, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "1.Drop_PHKA2_PLEKHA1.Cox.", run_id, ".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + EGFR + CP, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "Final.9Genes.Cox.", run_id, ".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# # manually test best model for V3 ------------------------------------------------
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ stage.numeric + EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + sex.ismale + EGFR + CP + PLEKHA1 + PHKA2 + TGFA, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "0.Drop_PLIN21.Cox.", run_id, ".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ stage.numeric + EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + sex.ismale + EGFR + CP + PLEKHA1 + PHKA2, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "1.Drop_PLIN21_TGFA.Cox.", run_id, ".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + EGFR + CP, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "Final.9Genes.Cox.", run_id, ".txt")
# sink(file2write)
# fit_efs_sum
# sink()
# 
# fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ stage.numeric + sex.ismale, data = testdata_df)
# fit_efs_sum <- summary(fit_efs)
# file2write <- paste0(dir_out, "Stage_Sex.Cox.", run_id, ".txt")
# sink(file2write)
# fit_efs_sum
# sink()

# plot --------------------------------------------------------------------
fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ stage.numeric + sex.ismale + EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + EGFR + CP, data = testdata_df)
file2write <- paste0(dir_out, "Hazard_Ratio.MultiGeneModel.pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
p <- ggforest(model = fit_efs, data = testdata_df)
print(p)
dev.off()

fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ stage.numeric + sex.ismale + EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + EGFR, data = testdata_df)
file2write <- paste0(dir_out, "Hazard_Ratio.MultiGeneModel.noCP.pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
p <- ggforest(model = fit_efs, data = testdata_df)
print(p)
dev.off()

fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ CP, data = testdata_df)
file2write <- paste0(dir_out, "Hazard_Ratio.CPonly.pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
p <- ggforest(model = fit_efs, data = testdata_df)
print(p)
dev.off()

fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~  sex.ismale + EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + EGFR + CP, data = testdata_df)
file2write <- paste0(dir_out, "Hazard_Ratio.MultiGeneModel.noStage.pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
p <- ggforest(model = fit_efs, data = testdata_df)
print(p)
dev.off()



fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~  EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + EGFR + CP, data = testdata_df)
file2write <- paste0(dir_out, "Hazard_Ratio.MultiGeneModel.ExpressionOnly.pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
p <- ggforest(model = fit_efs, data = testdata_df)
print(p)
dev.off()

fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ stage.numeric + sex.ismale, data = testdata_df)
file2write <- paste0(dir_out, "Hazard_Ratio.MultiGeneModel.Stage_Sex.pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
p <- ggforest(model = fit_efs, data = testdata_df)
print(p)
dev.off()

