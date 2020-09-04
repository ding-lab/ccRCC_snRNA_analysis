# Yige Wu @WashU Aug 2020
## 2020-09-03 removed the silent mutations

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
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input 10xmapping result
snRNA_mutation_df <- fread("./Resources/Analysis_Results/mutation/unite_10xmapping/20200303.v1/10XMapping.20200303.v1.tsv", data.table = F)

# map aliquot to case -----------------------------------------------------
snRNA_mutation_df$Case <- mapvalues(x = snRNA_mutation_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))

# filter to variant allele only --------------------------------------------
snRNA_mutation_filtered_df <- snRNA_mutation_df %>%
  filter(Case %in% c("C3L-00079", "C3L-00088", "C3N-01200")) %>%
  filter(allele_type == "Var") %>%
  select(Case, gene_symbol, mutation) %>%
  unique()

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Mutations_Mapped_to_snRNA.3cases.", "tsv")
write.table(x = snRNA_mutation_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
