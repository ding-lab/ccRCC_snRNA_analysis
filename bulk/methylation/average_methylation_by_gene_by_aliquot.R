# Yige Wu @WashU March 2020
## running on local
## for averaging the methylation beta values per gene per aliquot

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
## input the methylation by probe data
methylation_by_probe_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Bulk_Processed_Data/Methylation/CCRCC.TumorBeta.txt", data.table = F)

# transform and average ---------------------------------------------------
methylation_by_probe_df <- as.data.frame(methylation_by_probe_df)
table(colnames(methylation_by_probe_df))
## remove duplicate column names
methylation_by_probe_df <- methylation_by_probe_df[,!duplicated(colnames(methylation_by_probe_df))]
## add a column for gene symbol
methylation_by_probe_df <- methylation_by_probe_df %>%
  dplyr::mutate(gene_symbol = str_split_fixed(string = V1, pattern = "@", n = 3)[,3])
methylation_by_probe_df %>%
  select(gene_symbol) %>%
  table()
## get case ids
case_ids <- colnames(methylation_by_probe_df)[!(colnames(methylation_by_probe_df) %in% c("V1", "gene_symbol"))]
## group by gene and summarise by averaging
methylation_by_gene_df <- aggregate(x = methylation_by_probe_df[,case_ids], list(methylation_by_probe_df$gene_symbol), mean)
## rename column
methylation_by_gene_df <- methylation_by_gene_df %>%
  rename(gene_symbol = Group.1)
# write table -------------------------------------------------------------
file2write <- paste0(dir_out, "CCRCC.TumorBeta.", "by_gene.", run_id, ".tsv")
write.table(x = methylation_by_gene_df, file = file2write, sep = "\t", quote = F, row.names = F)
