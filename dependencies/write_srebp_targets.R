# Yige Wu @WashU Feb 2020
## for writing a table for NRF2 targets

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
## input TF integrations
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/PPI/TF_interactions.txt", data.table = F)

# filter ------------------------------------
srebp_tf_tab <- tf_tab %>%
  filter(source_genesymbol %in% c("SREBF1", "SREBF2")) %>%
  select(source_genesymbol, target_genesymbol)

# add some canonical targets from literature ------------------------------

# write table -------------------------------------------------------------
write.table(x = srebp_tf_tab, file = paste0(dir_out, "SREBP_Target_Genes.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
