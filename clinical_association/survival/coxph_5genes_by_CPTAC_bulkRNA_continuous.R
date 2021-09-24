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
library(RcppAlgos)
library(doParallel)
## set run id
version_tmp <- 2
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
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210824.v1/ccRCC_markers.Surface.20210824.v1.tsv")
genes_process <- genes_process_df$Gene
genes_process <- genes_process[!(genes_process %in% c("DPP6", "CPNE8", "EFNA5", "MGLL", "SPIRE1", "SPIRE1", "PLCB1", "OSMR", "SORBS1", "ANO6", "EPB41", "PAM", "RHEX"))]
gene1234_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/clinical_association/survival/coxph_4genes_by_CPTAC_bulkRNA_continuous/20210920.v2/Cox_4genes.bulkRNA.continuous.Filtered.20210920.v2.tsv")

result_list <- list()
# for (gene_test in "CP") {
for (i in 1:nrow(gene1234_process_df)) {
  result_list[[i]] <- list()
  gene1_tmp <- gene1234_process_df$genesymbol1[i]
  gene2_tmp <- gene1234_process_df$genesymbol2[i]
  gene3_tmp <- gene1234_process_df$genesymbol3[i]
  gene4_tmp <- gene1234_process_df$genesymbol4[i]
  
  gene5_process <- genes_process[!(genes_process %in% c(gene1_tmp, gene2_tmp, gene3_tmp, gene4_tmp))]
  for (gene5_tmp in gene5_process) {
    genes_process_tmp <- c(gene1_tmp, gene2_tmp, gene3_tmp, gene4_tmp, gene5_tmp)
    exp_test_wide_df <- exp_data_df[exp_df$gene_name %in% genes_process_tmp,]
    testdata_df <- data.frame(t(exp_test_wide_df))
    
    genes_process_colnames <- exp_df$gene_name[exp_df$gene_name %in% genes_process_tmp]
    colnames(testdata_df)[genes_process_colnames == gene1_tmp] <- "Expression.gene1"
    colnames(testdata_df)[genes_process_colnames == gene2_tmp] <- "Expression.gene2"
    colnames(testdata_df)[genes_process_colnames == gene3_tmp] <- "Expression.gene3"
    colnames(testdata_df)[genes_process_colnames == gene4_tmp] <- "Expression.gene4"
    colnames(testdata_df)[genes_process_colnames == gene5_tmp] <- "Expression.gene5"
    
    testdata_df$CASE_ID <- colnames(exp_test_wide_df)
    testdata_df <- merge(x = testdata_df, y = survival_df, by = c("CASE_ID"), all.x = T)
    
    testdata_df <- testdata_df %>%
      filter(CASE_ID != "C3L-00359") %>%
      mutate(EFS_censor = (with_new_event == "Tumor Free")) %>%
      mutate(EFS = (survival_time + 9)/365) %>%
      arrange(CASE_ID, desc(Expression.gene1))
    testdata_df <- testdata_df[!duplicated(testdata_df$CASE_ID),]
    
    ## EFS_censor == 0 with event; == 1 without event
    ## test
    testdata_comp_df <- testdata_df %>% 
      filter(!is.na(EFS_censor) & !is.na(EFS))
    
    fit_efs <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ Expression.gene1 + Expression.gene2 + Expression.gene3 + Expression.gene4 + Expression.gene5, data = testdata_comp_df)
    fit_efs_sum <- summary(fit_efs)
    result_list[[i]][[gene5_tmp]] <- c(genes_process_tmp,
                  fit_efs_sum$waldtest["pvalue"], fit_efs_sum$logtest["pvalue"], fit_efs_sum$sctest["pvalue"], 
                  fit_efs_sum$concordance[1],fit_efs_sum$concordance[2],
                  fit_efs_sum$coefficients[1,1], fit_efs_sum$coefficients[2,1], fit_efs_sum$coefficients[3,1], fit_efs_sum$coefficients[4,1], fit_efs_sum$coefficients[5,1])
  }
}

cox_result_df <- do.call(cbind.data.frame, result_list)
cox_result_df <- data.frame(t(cox_result_df))

colnames(cox_result_df) <- c("genesymbol1", "genesymbol2", "genesymbol3", "genesymbol4", "genesymbol5",
                             "pvalue.wald", "pvalue.lr", "pvalue.logrank", 
                             "concordance", "se.concordance",
                             "coef.gene1", "coef.gene2", "coef.gene3", "coef.gene4", "coef.gene5")
cox_result_df$concordance <- as.numeric(cox_result_df$concordance)
cox_result_df$se.concordance <- as.numeric(cox_result_df$se.concordance)

cox_result_df <- cox_result_df %>%
  mutate(CI.low.concordance = concordance - 1.96*se.concordance)

cox_result_filtered_df <- cox_result_df %>%
  filter(pvalue.logrank < 0.05) %>%
  arrange(desc(concordance))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cox_5genes.bulkRNA.continuous.", run_id, ".tsv")
write.table(x = cox_result_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Cox_5genes.bulkRNA.continuous.Filtered.", run_id, ".tsv")
write.table(x = cox_result_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
