# Yige Wu @WashU Dec 2020
## 3 categories
### partial 3p  loss: require for each CNV, the cluster with the most %cells with CNV sould be over 10% while the cluster with the least %cells with CNVs should be M10%

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
cnvtype_3ploss_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/assign_3ploss_status_by_tumor/20201207.v1/3pLoss_Type_Assignment_Per_Tumor.20201207.v1.tsv")
cnvtype_5qgain_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/assign_5qgain_status_by_tumor/20201207.v1/5qGain_Type_Assignment_Per_Tumor.20201207.v1.tsv")

# combine -----------------------------------------------------------------
cnvtype_merged_df <- merge(x = cnvtype_3ploss_df, y = cnvtype_5qgain_df, by = c("aliquot.wu"), all = T)
cnvtype_count_df <- cnvtype_merged_df %>%
  select(Group_3ploss, Group_5qgain) %>%
  table()
cnvtype_count_df

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CNVType_3pLoss_5qGain_Per_Tumor.", run_id, ".tsv")
write.table(x = cnvtype_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "CNVType_Count_3pLoss_5qGain.", run_id, ".tsv")
write.table(x = cnvtype_count_df, file = file2write, quote = F, sep = "\t", row.names = T)
