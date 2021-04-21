# Yige Wu @WashU Oct 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv", data.table = F)
## input the spearman pairwise correlation result
pearson_coef.tumorcellvariable_genes.df <- fread(input = "./Resources/Analysis_Results/pairwise_correlation/calculate_tumor_manualsubcluster_pairwise_correlation_tumorcellvariable_genes/20201201.v1/avg_exp_by_tumorsubluster.tumorcellvaraible_genes.pearson_coef20201201.v1.tsv", data.table = F)

# get ids -----------------------------------------------------------------
ids_tumorsubcluster <- pearson_coef.tumorcellvariable_genes.df$V1
names_tumorsubclusters <- gsub(x = ids_tumorsubcluster, pattern = "\\.", replacement = "-")
names_tumorsubclusters
## get reable sample id
ids_aliquot_wu <- str_split_fixed(string = names_tumorsubclusters, pattern = "_", n = 2)[,1]
case_ids <- str_split_fixed(string = names_tumorsubclusters, pattern = "-T", n = 2)[,1]

# process pairwise correlation -----------------------------------------------
coef_long_df <- melt(data = pearson_coef.tumorcellvariable_genes.df)
coef_long_filtered_df <- coef_long_df %>%
  filter(value < 1) %>%
  rename(id_tumorsubcluster1 = V1) %>%
  rename(id_tumorsubcluster2 = variable)
coef_long_filtered_df$id_case1 <- mapvalues(x = coef_long_filtered_df$id_tumorsubcluster1, from = ids_tumorsubcluster, to = case_ids)
coef_long_filtered_df$id_case2 <- mapvalues(x = coef_long_filtered_df$id_tumorsubcluster2, from = ids_tumorsubcluster, to = case_ids)
## filter to only correlations within the same patient
coef_long_filtered_df2 <- coef_long_filtered_df %>%
  filter(id_case1 == id_case2)
## summarize the average correlation per case
coef_summary_df <- coef_long_filtered_df2 %>%
  group_by(id_case1) %>%
  summarise(coef_mean = mean(value)) %>%
  mutate(ith_score = (1-coef_mean)) %>%
  arrange(desc(ith_score))

# write outupt ------------------------------------------------------------
file2write <- paste0(dir_out, "Mean_Pairwise_Corrleation_TumorSubclusters_By_Case.", run_id, ".tsv")
write.table(x = coef_summary_df, file = file2write, sep = "\t", row.names = F, quote = F)

