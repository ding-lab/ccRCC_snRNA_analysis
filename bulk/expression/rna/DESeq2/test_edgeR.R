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
library(data.table)
library(dplyr)
library(plyr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input mutation table
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210310.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210310.v1.tsv")
## input sample info
sampleinfo_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/gdc_sample_sheet.2021-03-11.tsv")
## input htseq-count file info
# dir_htseq <- "./Resources/Bulk_Processed_Data/mRNA/HTSeq-count.gdc_download_20210311_143816.655582/"
dir_htseq <- "./Resources/Bulk_Processed_Data/mRNA/HTSeq-count-gzfiles/"

# preprocess --------------------------------------------------------------
sampleinfo_filtered_df <- sampleinfo_df %>%
  mutate(Case = str_split_fixed(string = `Case ID`, pattern = ",", n = 2)[,1]) %>%
  filter(Case %in% mut_df$Case)
sampleinfo_filtered_tumor_df <- sampleinfo_filtered_df %>%
  filter(`Sample Type` == "Primary Tumor")
sampleinfo_filtered_nat_df <- sampleinfo_filtered_df %>%
  filter(`Sample Type` == "Solid Tissue Normal") %>%
  mutate(condition = "NAT")
## add mutation status
sampleinfo_filtered_tumor_df$condition <- mapvalues(x = sampleinfo_filtered_tumor_df$Case, from = mut_df$Case, to = as.vector(mut_df$mutation_category_sim))
## combine
sampleinfo_filtered_df <- rbind(sampleinfo_filtered_tumor_df, sampleinfo_filtered_nat_df)


# Creating a DGEList object -----------------------------------------------
# dgelist_obj <- readDGE(files = "0a1a14f8-3e7a-44ae-a7ef-c0ae33a66975/8a95b044-5ec9-45b6-8e33-ed2a12d4b547.htseq_counts.txt.gz", path = dir_htseq) ## this no work
# dgelist_obj <- readDGE(files = "0a1a14f8-3e7a-44ae-a7ef-c0ae33a66975/8a95b044-5ec9-45b6-8e33-ed2a12d4b547.htseq_counts.txt", path = dir_htseq) ## this works
# dgelist_obj_test <- readDGE(files = "8a95b044-5ec9-45b6-8e33-ed2a12d4b547.htseq_counts.txt.gz", path = paste0(dir_htseq, "0a1a14f8-3e7a-44ae-a7ef-c0ae33a66975/"), labels = "test") ## this works
dgelist_obj <- readDGE(files = sampleinfo_filtered_df$`File Name`, path = paste0(dir_htseq), labels = sampleinfo_filtered_df$`Sample ID`) ## this works
dgelist_obj$counts %>% nrow() # [1] 60487

# Filtering ---------------------------------------------------------------
countsPerMillion <- cpm(dgelist_obj)
summary(countsPerMillion)
countCheck <- countsPerMillion > 1
rm(countsPerMillion)
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
rm(countCheck)
dgelist_obj <- dgelist_obj[keep,]
dim(dgelist_obj) ## [1] 22888   257
summary(cpm(dgelist_obj))

# Normalisation -----------------------------------------------------------
dgelist_obj <- calcNormFactors(dgelist_obj, method="TMM")

# Data Exploration --------------------------------------------------------
file2write <- paste0(dir_out, "MDS.png")
png(file2write, width = 2000, height = 2000, res = 150)
plotMDS(dgelist_obj)
dev.off()

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
fit <- glmFit(dgelist_obj, designMat)
View(fit$design)
lrt_bap1vsnat <- glmLRT(fit, coef=2)
result_bap1vsnat <- topTags(lrt_bap1vsnat, n = 15000)
View(result_bap1vsnat$table)
deGenes <- decideTestsDGE(lrt_bap1vsnat, p=0.001)
deGenes <- rownames(lrt_bap1vsnat)[as.logical(deGenes)]
plotSmear(lrt_bap1vsnat, de.tags=deGenes)
abline(h=c(-1, 1), col=2)

