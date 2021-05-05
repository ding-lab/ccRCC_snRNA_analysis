# Yige Wu @WashU May 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/load_data.R")

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input bulk mutation result
maf_discovery_df1 <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Somatic_mutation/CPTAC_ccRCC_discovery_somatic_mutation.maf_v1.0.tsv")
maf_discovery_df2 <- loadMaf()
## input bulk mutation result for the ITH cohort
maf_ith_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_ITH/Somatic_mutation/CPTAC_ccRCC_ITH_somatic_mutation.maf_v1.0.tsv")
## input id meta data
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv", data.table = F)

# format somatic mutation data for the discovery data ----------------------------------------------------
## set case ids to process
case_ids <- unique(id_metadata_df$Case[id_metadata_df$snRNA_available])
sample_ids <- unique(id_metadata_df$Sample[id_metadata_df$snRNA_available])
## filter
colnames_int <- intersect(colnames(maf_discovery_df1), colnames(maf_discovery_df2))
maf_discovery_df <- rbind(maf_discovery_df1[,colnames_int], maf_discovery_df2[,colnames_int]) %>% unique()
maf_discovery_filtered_df <- maf_discovery_df %>%
  mutate(case_id = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", n = 2)[,1]) %>%
  filter(case_id %in% case_ids)
maf_discovery_filtered_df$sample_id <- mapvalues(x = maf_discovery_filtered_df$case_id, from = subset(x = id_metadata_df, Is_discovery_set == T & Sample_Type == "Tumor")$Case, to = as.vector(subset(x = id_metadata_df, Is_discovery_set == T & Sample_Type == "Tumor")$Sample))
maf_ith_filtered_df <- maf_ith_df %>%
  mutate(sample_id = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", n = 2)[,1]) %>%
  filter(sample_id %in% sample_ids)
colnames_int <- intersect(colnames(maf_ith_filtered_df), colnames(maf_discovery_filtered_df))
maf_merged_df <- rbind(maf_discovery_filtered_df[, colnames_int], 
                       maf_ith_filtered_df[, colnames_int]) %>% unique()

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "snRNA_samples.somatic_mutation.maf.v1.0.tsv")
write.table(x = maf_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
