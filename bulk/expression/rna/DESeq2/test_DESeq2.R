# Yige Wu @WashU Mar 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
# BiocManager::install("DESeq2")
library(DESeq2)
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
dir_htseq <- "./Resources/Bulk_Processed_Data/mRNA/HTSeq-count.gdc_download_20210311_143816.655582/"

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
## transform
sampleTable <- sampleinfo_filtered_df %>%
  mutate(sampleName = paste0(`File ID`, "/", `File Name`)) %>%
  mutate(fileName = sampleName) %>%
  select(sampleName, fileName, condition)

# build DESeqDataSet object -----------------------------------------------
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

