# Yige Wu @WashU March 2020
## for calculating how many manually grouped tumor subclusters each tumor has

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
## input manual subcluster grouping
manual_tumorsubcluster_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/Tumor_Subcluster/ccRCC_snRNA_Downstream_Processing - Individual.TumorCluster2Cell_Type.20200223.v1.tsv", data.table = F)

# count the number of manual subclusters per sample ----------------------
num_manual_tumorsubcluster_df <- manual_tumorsubcluster_df %>%
  filter(Cluster_Manual != "") %>%
  select(Aliquot, Cluster_Manual) %>%
  unique() %>%
  select(Aliquot) %>%
  table() %>%
  as.data.frame() %>%
  arrange(-Freq) %>%
  rename(Aliquot = '.') %>%
  rename(Number_Manual_TumorClusters = Freq)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "number_of_manual_tumorsubcluster_per_sample.", run_id, ".tsv")
write.table(x = num_manual_tumorsubcluster_df, file = file2write, quote = F, sep = "\t", row.names = F)

