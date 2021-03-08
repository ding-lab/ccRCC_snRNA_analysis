# Yige Wu @WashU March 2020
## merge the tumor content estimated from snRNA data with bulk RNA data

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the tumor content estimate from the snRNA data
snRNA_tumorcontent_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/other/calculate_tumor_cell_fraction/20200218.v1/Tumor_Cell_Fraction.20200218.v1.tsv", data.table = F)
## only keep the tumor samples
snRNA_tumorcontent_df <- snRNA_tumorcontent_df %>%
  filter(Sample_Type == "Tumor")
## input the tumor content estimate from the bulk RNA data (ESTIMATE)
estimate_df <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S7.xlsx", sheet = "ESTIMATE scores")
## input meta data
id_metadata_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)

# extract the bulk RNA-based tumor content estimate -----------------------
bulkRNA_tumor_content <- as.numeric(as.vector(as.data.frame(estimate_df[7,2:176])))
bulkRNA_tumor_content_df <- data.frame(ESTIMATE_TumorPurity_RNA = bulkRNA_tumor_content, Specimen.ID.bulk = unlist(estimate_df[3,2:176]))
bulkRNA_tumor_content_df <- bulkRNA_tumor_content_df %>%
  filter(Specimen.ID.bulk %in% id_metadata_df$Aliquot.bulk)
## map bulk ID to aliquot ids for snRNA-Seq
bulkRNA_tumor_content_df$Specimen.ID.snRNA <- mapvalues(bulkRNA_tumor_content_df$Specimen.ID.bulk, from = id_metadata_df$Aliquot.bulk, to = as.vector(id_metadata_df$Aliquot.snRNA))
## filter by snRNA aliquot id: only keep those with available snRNA data
bulkRNA_tumor_content_df <- bulkRNA_tumor_content_df %>%
  filter(Specimen.ID.snRNA %in% snRNA_tumorcontent_df$Aliquot)

# merge ESTIMATE tumor content into sn-based tumor content ------------------------------
## only keep the tumor samples
merged_tumorcontent_df <- snRNA_tumorcontent_df
merged_tumorcontent_df$ESTIMATE_TumorPurity_RNA <- mapvalues(merged_tumorcontent_df$Aliquot, from = bulkRNA_tumor_content_df$Specimen.ID.snRNA, to = as.vector(bulkRNA_tumor_content_df$ESTIMATE_TumorPurity_RNA))
## filter out samples without bulk RNA-based deconvolution
merged_tumorcontent_df$ESTIMATE_TumorPurity_RNA[merged_tumorcontent_df$ESTIMATE_TumorPurity_RNA == merged_tumorcontent_df$Aliquot] <- NA

# write table -------------------------------------------------------------
file2write <- paste0(dir_out, "Perc_Tumor_Content_from_snRNA_and_bulkRNA.", run_id, ".tsv")
write.table(x = merged_tumorcontent_df, file = file2write, quote = F, row.names = F, sep = "\t")

