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
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210824.v1/ccRCC_markers.Surface.20210824.v1.tsv")
genes_process <- genes_process_df$Gene
genes_process <- genes_process[!(genes_process %in% c("DPP6", "CPNE8", "EFNA5", "MGLL", "SPIRE1", "SPIRE1", "PLCB1", "OSMR", "SORBS1", "ANO6", "EPB41", "PAM", "RHEX"))]

result_list <- list()
# for (gene_test in "CP") {
for (gene_test in genes_process) {
  # make combined data and test ------------------------------------------------------
  ## filter specific protein data
  exp_test_wide_df <- exp_data_df[exp_df$V1 == gene_test,]
  testdata_df <- data.frame(CASE_ID = exp_colnames_filtered_df$CASE_ID, Expression = unlist(exp_test_wide_df))
  testdata_df <- merge(x = testdata_df, y = survival_df, by = c("CASE_ID"), all.x = T)
  cutoff_exp_low <- quantile(x = testdata_df$Expression, probs = 0.35, na.rm = T); cutoff_exp_low
  cutoff_exp_high <- quantile(x = testdata_df$Expression, probs = 0.65, na.rm = T); cutoff_exp_high
  
  testdata_df <- testdata_df %>%
    filter(CASE_ID != "C3L-00359") %>%
    mutate(Expression_group = ifelse(Expression < cutoff_exp_low, "Low", ifelse(Expression > cutoff_exp_high, "High", "Medium"))) %>%
    mutate(EFS_censor = (with_new_event == "Tumor Free")) %>%
    mutate(EFS = (survival_time + 9)/365) %>%
    arrange(CASE_ID, desc(Expression))
  testdata_df <- testdata_df[!duplicated(testdata_df$CASE_ID),]
  
  ## EFS_censor == 0 with event; == 1 without event
  ## test
  testdata_comp_df <- testdata_df %>%
    filter(!is.na(EFS_censor) & !is.na(EFS) & !is.na(Expression_group))
  testdata_comp_df$Expression_group <- factor(x = testdata_comp_df$Expression_group, levels = c("High", "Medium", "Low"))
  
  fit_efs <- coxph(formula = Surv(EFS, EFS_censor == 0) ~ Expression, data = testdata_comp_df)
  fit_efs_sum <- summary(fit_efs)
  result_list[[gene_test]] <- c(gene_test,
                                fit_efs_sum$waldtest["pvalue"], fit_efs_sum$logtest["pvalue"], fit_efs_sum$sctest["pvalue"], 
                                fit_efs_sum$concordance[1],
                                fit_efs_sum$coefficients[1, "coef"])
}
cox_result_df <- do.call(rbind.data.frame, result_list)
colnames(cox_result_df) <- c("genesymbol", "pvalue.wald", "pvalue.lr", "pvalue.logrank", "concordance", "coef")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cox_1gene.continuous_variable.", run_id, ".tsv")
write.table(x = cox_result_df, file = file2write, quote = F, sep = "\t", row.names = F)
