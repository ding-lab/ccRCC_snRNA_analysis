# Yige Wu @WashU Jun 2021
## test for overall survival

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
library(RcppAlgos)
library(doParallel)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input RNA-seq RSEM normalized data, needs to be log-transformed
### more about the input data: 
#### https://www.cbioportal.org/faq
#### https://www.biostars.org/p/106127/
exp_df <- fread(data.table = F, input = "~/Box/Ding_Lab/Databases/cBioPortal/CCRCC/kirc_tcga_pub/data_RNA_Seq_v2_expression_median.txt")
## input survival ddata
survival_df <- fread(data.table = F, input = "~/Box/Ding_Lab/Databases/cBioPortal/CCRCC/kirc_tcga_pan_can_atlas_2018_clinical_data.tsv")

# specify gene to test ----------------------------------------------------
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210824.v1/ccRCC_markers.Surface.20210824.v1.tsv")
genes_process <- genes_process_df$Gene
genes_process <- genes_process[!(genes_process %in% c("DPP6", "CPNE8", "EFNA5", "MGLL", "SPIRE1", "SPIRE1", "PLCB1", "OSMR", "SORBS1", "ANO6", "EPB41", "PAM", "RHEX"))]

# preprocess ------------------------------------------------------
sampleids <- colnames(exp_df)
sampleids <- sampleids[3:length(sampleids)]
exp_data_df <- as.matrix(exp_df[,sampleids])
exp_data_df <- log2(exp_data_df+1)
exp_data_df <- as.data.frame(exp_data_df)
dim(exp_data_df)
dim(exp_df)
##
exp_test_wide_df <- exp_data_df[exp_df$Hugo_Symbol %in% genes_process,]
testdata_df <- data.frame(t(exp_test_wide_df))
colnames(testdata_df) <- exp_df$Hugo_Symbol[exp_df$Hugo_Symbol %in% genes_process]
testdata_df$CASE_ID <- gsub(x = colnames(exp_test_wide_df), pattern = "\\-01", replacement = "")
testdata_df <- merge(x = testdata_df, y = survival_df %>%
                       rename(CASE_ID = 'Patient ID') %>%
                       rename(tumor_stage_pathological = 'Neoplasm Disease Stage American Joint Committee on Cancer Code') %>%
                       mutate(tumor_stage_pathological = gsub(x = tumor_stage_pathological, pattern = "Stage ", replacement = "")) %>%
                       mutate(stage.numeric = ifelse(tumor_stage_pathological == "I", 1,
                                                     ifelse(tumor_stage_pathological == "II", 2,
                                                            ifelse(tumor_stage_pathological == "III", 3, 4)))) %>%
                       mutate(sex.ismale = as.numeric(Sex == "Male")) %>%
                       rename(EFS_month = `Progress Free Survival (Months)`) %>%
                       mutate(EFS = EFS_month/12) %>%
                       rename(EFS_censor = `Progression Free Status`) %>%
                       mutate(EFS_censor = as.numeric(str_split_fixed(string = EFS_censor, pattern = ":", n = 2)[,1])) %>%
                       select(CASE_ID, stage.numeric, sex.ismale, EFS, EFS_censor), by = "CASE_ID", all.x = T)
testdata_df <- testdata_df %>%
  filter(!is.na(EFS_censor) & !is.na(EFS))


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

fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ EPHA6 + ABCC3 + FTO + COL23A1 + CA9 + SEMA6A + NDRG1 + EGFR + CP, data = testdata_df)
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

