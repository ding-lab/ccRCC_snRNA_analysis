# Yige Wu @WashU Mar 2021
## reference: https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
if (!requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR")
}
library(edgeR)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input mutation table
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210310.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210310.v1.tsv")
## input meta data
ccrcc_df <- readxl::read_excel(path = "./Resources/Bulk_Processed_Data/CPTAC3-ccRCC-SupplementaryTables_Final/Table S1.xlsx", sheet = "ccrcc_clinical_characteristics")
## input sample info
sampleinfo_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/gdc_sample_sheet.2021-03-11.tsv")
## input the DEGList object
dgelist_obj <- readRDS(file = './Resources/Analysis_Results/bulk/expression/edgeR/create_DGEList_object/create_DEGList_CPTAC_ccRCC_Discovery_cases/20210312.v1/CPTAC_Discovery_ccRCC_Cases.DEGList.RDS')
## input the gene id mapping
gene_mapping_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/map_to_genesymbols/map_CPTAC_Discovery_ensemblgeneids_to_genesymbols/20210312.v1/ensembl_gene_id.mapping.nodup.20210312.v1.tsv")

# preprocess --------------------------------------------------------------
sampleinfo_filtered_df <- sampleinfo_df %>%
  mutate(Case = str_split_fixed(string = `Case ID`, pattern = ",", n = 2)[,1]) %>%
  filter(Case %in% ccrcc_df$Case_ID[ccrcc_df$Histologic_Type == "Clear cell renal cell carcinoma"])
sampleinfo_filtered_tumor_df <- sampleinfo_filtered_df %>%
  filter(`Sample Type` == "Primary Tumor")
sampleinfo_filtered_nat_df <- sampleinfo_filtered_df %>%
  filter(`Sample Type` == "Solid Tissue Normal") %>%
  mutate(condition = "NAT")
## add mutation status
sampleinfo_filtered_tumor_df$condition <- mapvalues(x = sampleinfo_filtered_tumor_df$Case, from = mut_df$Case, to = as.vector(mut_df$mutation_category_sim))
sampleinfo_filtered_tumor_df$condition[as.vector(sampleinfo_filtered_tumor_df$condition) == as.vector(sampleinfo_filtered_tumor_df$Case)] <- "Non-mutants"
## combine
sampleinfo_filtered_df <- rbind(sampleinfo_filtered_tumor_df, sampleinfo_filtered_nat_df)

# Setting up the Model ----------------------------------------------------
sampleGroup <- mapvalues(x = colnames(dgelist_obj), from = sampleinfo_filtered_df$`Sample ID`, to = as.vector(sampleinfo_filtered_df$condition))
sampleGroup <- factor(x = sampleGroup, levels = c("NAT", "BAP1 mutated", "PBRM1 mutated", "Non-mutants", "Both mutated"))
designMat <- model.matrix(~sampleGroup)
View(designMat)

#  Estimating Dispersions -------------------------------------------------
dgelist_obj <- estimateGLMCommonDisp(dgelist_obj, design=designMat)
dgelist_obj <- estimateGLMTrendedDisp(dgelist_obj, design=designMat)
dgelist_obj <- estimateGLMTagwiseDisp(dgelist_obj, design=designMat)
file2write <- paste0(dir_out, "BCV.pdf")
pdf(file2write, width = 15, height = 15)
plotBCV(dgelist_obj)
dev.off()

# 2.8 Differential Expression ---------------------------------------------
fit <- glmQLFit(dgelist_obj, designMat)

## BAP1 vs Non-mutants + PBRM1
qlf <- glmQLFTest(fit, contrast = c(0, 1, -1, -1, 0))
result_df <- topTags(object = qlf, n = nrow(qlf$table))$table
result_df$ensembl_gene_id_version.deg <- rownames(result_df)

# map gene symbol ---------------------------------------------------------
result_df <- result_df %>%
  mutate(ensembl_gene_id = str_split_fixed(string = ensembl_gene_id_version.deg, pattern = "\\.", n = 2)[,1])

result_df <- merge(x = result_df, y = gene_mapping_df, by = c("ensembl_gene_id"), all.x = T)
result_df <- result_df %>%
  dplyr::filter(!(ensembl_gene_id_version.deg %in% c("__alignment_not_unique", "__ambiguous", "__no_feature")))

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_Mutated_vs_PBRM1_NonMutants.glmQLFTest.OutputTables.tsv")
write.table(x = result_df, file = file2write, sep = "\t", row.names = F, quote = F)



