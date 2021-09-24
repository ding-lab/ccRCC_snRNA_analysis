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
library(TCGAbiolinks)
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
survival_df <- fread(data.table = F, input = "~/Box/Ding_Lab/Databases/cBioPortal/CCRCC/kirc_tcga_pub/data_clinical_patient.txt", skip = 4)
clinical_df <- GDCquery_clinic(project = "TCGA-KIRC", type = "clinical", save.csv = FALSE)
survival_df <- fread(data.table = F, input = "~/Box/Ding_Lab/Databases/cBioPortal/CCRCC/kirc_tcga_pan_can_atlas_2018_clinical_data.tsv")


# preprocess ------------------------------------------------------
sampleids <- colnames(exp_df)
sampleids <- sampleids[3:length(sampleids)]
exp_data_df <- as.matrix(exp_df[,sampleids])
exp_data_df <- log2(exp_data_df+1)
exp_data_df <- as.data.frame(exp_data_df)
dim(exp_data_df)
dim(exp_df)

# specify gene to test ----------------------------------------------------
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210824.v1/ccRCC_markers.Surface.20210824.v1.tsv")
genes_process <- genes_process_df$Gene
genes_process <- genes_process[!(genes_process %in% c("DPP6", "CPNE8", "EFNA5", "MGLL", "SPIRE1", "SPIRE1", "PLCB1", "OSMR", "SORBS1", "ANO6", "EPB41", "PAM", "RHEX"))]
gene123_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/clinical_association/survival/coxph_3genes_by_CPTAC_bulkRNA_continuous/20210920.v1/Cox_3genes.bulkRNA.continuous.Filtered.20210920.v1.tsv")

result_list <- list()
# for (gene_test in "CP") {
for (i in 1:nrow(gene123_process_df)) {
  result_list[[i]] <- list()
  gene1_tmp <- gene123_process_df$genesymbol1[i]
  gene2_tmp <- gene123_process_df$genesymbol2[i]
  gene3_tmp <- gene123_process_df$genesymbol3[i]
  
  gene4_process <- genes_process[!(genes_process %in% c(gene1_tmp, gene2_tmp, gene3_tmp))]
  for (gene4_tmp in gene4_process) {
    genes_process_tmp <- c(gene1_tmp, gene2_tmp, gene3_tmp, gene4_tmp)
    exp_test_wide_df <- exp_data_df[exp_df$Hugo_Symbol %in% genes_process_tmp,]
    testdata_df <- data.frame(t(exp_test_wide_df))
    
    genes_process_colnames <- exp_df$Hugo_Symbol[exp_df$Hugo_Symbol %in% genes_process_tmp]
    colnames(testdata_df)[genes_process_colnames == gene1_tmp] <- "Expression.gene1"
    colnames(testdata_df)[genes_process_colnames == gene2_tmp] <- "Expression.gene2"
    colnames(testdata_df)[genes_process_colnames == gene3_tmp] <- "Expression.gene3"
    colnames(testdata_df)[genes_process_colnames == gene4_tmp] <- "Expression.gene4"
    
    testdata_df$CASE_ID <- gsub(x = colnames(exp_test_wide_df), pattern = "\\-01", replacement = "")
    testdata_df <- merge(x = testdata_df, y = survival_df, by.x = c("CASE_ID"), by.y = c("PATIENT_ID"), all.x = T)
    
    testdata_df <- testdata_df %>%
      mutate(EFS_censor = (OS_STATUS == "DECEASED")) %>%
      mutate(EFS = OS_MONTHS) %>%
      arrange(CASE_ID, desc(Expression.gene1))
    testdata_df <- testdata_df[!duplicated(testdata_df$CASE_ID),]
    
    ## EFS_censor == 0 with event; == 1 without event
    ## test
    testdata_comp_df <- testdata_df %>% 
      filter(!is.na(EFS_censor) & !is.na(EFS))
    
    fit_efs <- coxph(formula = Surv(EFS, EFS_censor) ~ Expression.gene1 + Expression.gene2 + Expression.gene3 + Expression.gene4, data = testdata_comp_df)
    fit_efs_sum <- summary(fit_efs)
    result_list[[i]][[gene4_tmp]] <- c(genes_process_tmp,
                  fit_efs_sum$waldtest["pvalue"], fit_efs_sum$logtest["pvalue"], fit_efs_sum$sctest["pvalue"], 
                  fit_efs_sum$concordance[1],fit_efs_sum$concordance[2],
                  fit_efs_sum$coefficients[1,1], fit_efs_sum$coefficients[2,1], fit_efs_sum$coefficients[3,1], fit_efs_sum$coefficients[4,1])
  }
}

cox_result_df <- do.call(cbind.data.frame, result_list)
cox_result_df <- data.frame(t(cox_result_df))

colnames(cox_result_df) <- c("genesymbol1", "genesymbol2", "genesymbol3", "genesymbol4",
                             "pvalue.wald", "pvalue.lr", "pvalue.logrank", 
                             "concordance", "se.concordance",
                             "coef.gene1", "coef.gene2", "coef.gene3", "coef.gene4")
cox_result_df$concordance <- as.numeric(cox_result_df$concordance)
cox_result_df$se.concordance <- as.numeric(cox_result_df$se.concordance)

cox_result_df <- cox_result_df %>%
  mutate(CI.low.concordance = concordance - 1.96*se.concordance)

cox_result_filtered_df <- cox_result_df %>%
  filter(pvalue.logrank < 0.05) %>%
  filter(coef.gene1 > 0 & coef.gene2 > 0 & coef.gene3 > 0)
  # filter(coef.gene1 > 0 & coef.gene2 > 0 & coef.gene3 > 0 & coef.gene4 > 0)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cox_4genes.bulkRNA.continuous.", run_id, ".tsv")
write.table(x = cox_result_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Cox_4genes.bulkRNA.continuous.Filtered.", run_id, ".tsv")
write.table(x = cox_result_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
