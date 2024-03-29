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
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltype13_bysample_katmai/20210701.v1/33_aliquot_merged.avgexp.SCT.data.bycelltype13_bysample.20210701.v1.tsv")
## input bulk meta data
metadata_df <- fread("./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input survival ddata
survival_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_cptac_discovery_ccRCC_survival_time/20210621.v1/CPTAC_Discovery_ccRCC_Survival_Time20210621.v1.tsv")
survival_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_cptac_discovery_ccRCC_survival_time/20220317.v1/CPTAC_Discovery_ccRCC_Survival_Time20220317.v1.tsv")

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
# namescolors_expgroup <- c("High", "Low", "Medium")
colors_expgroup <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")[c(9, 4, 6)]

# specify gene to test ----------------------------------------------------
os_pvalue_vec <- NULL
pfs_pvalue_vec <- NULL
fontsize_plot <- 14
gene_test <- "PCSK6"
# make combined data and test ------------------------------------------------------
## filter specific protein data
exp_test_wide_df <- exp_data_df[exp_df$V1 == gene_test,]
testdata_df <- data.frame(CASE_ID = exp_colnames_filtered_df$CASE_ID, Expression = unlist(exp_test_wide_df))
testdata_df <- merge(x = testdata_df, y = survival_df, by = c("CASE_ID"), all.x = T)
cutoff_exp_low <- quantile(x = testdata_df$Expression, probs = 0.35, na.rm = T); cutoff_exp_low
cutoff_exp_high <- quantile(x = testdata_df$Expression, probs = 0.65, na.rm = T); cutoff_exp_high

testdata_df <- testdata_df %>%
  filter(CASE_ID != "C3L-00359") %>%
  mutate(Expression_group = ifelse(Expression <= cutoff_exp_low, "Low", ifelse(Expression >= cutoff_exp_high, "High", "Medium"))) %>%
  arrange(CASE_ID, desc(Expression))
testdata_df <- testdata_df[!duplicated(testdata_df$CASE_ID),]

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
                  legend.title = paste0(gene_test, " expression\n(snRNA-seq)"),
                  legend.labs = c("High", "Low"),
                  legend = "top",
                  xlab = "Time (years)",
                  ylab = "Overall Survival",
                  palette = c("#800026", "#FEB24C"),
                  ggtheme = theme_survminer(base_size = fontsize_plot,
                                            base_family = "",
                                            font.main = c(fontsize_plot, "plain", "black"),
                                            font.submain = c(fontsize_plot, "plain", "black"),
                                            font.x = c(fontsize_plot, "plain", "black"),
                                            font.y = c(fontsize_plot, "plain", "black"),
                                            font.caption = c(fontsize_plot, "plain", "black"),
                                            font.tickslab = c(fontsize_plot, "plain", "black"),
                                            legend = c("top", "bottom", "left", "right", "none"),
                                            font.legend = c(fontsize_plot, "plain", "black")),
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

## write source data
source_data_df <- testdata_comp_df %>%
  select(surv_time, surv_status, Expression_group)
write.table(x = source_data_df, file = paste0("~/Desktop/SF2f.", gene_test, ".snRNA.SourceData.tsv"), quote = F, sep = "\t", row.names = F)

