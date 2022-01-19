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
## input NRF target genes from Hayes et al. 2009
nrf_tf_manual1 <- read_excel(path = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/PPI/NRF2_target_genes.xlsx", sheet = "Hayes et al. 2009")
## input NRF target genes from Chorley et al. 2012
nrf_tf_manual2 <- read_excel(path = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/PPI/NRF2_target_genes.xlsx", sheet = "Chorley et al. 2012")

# filter ------------------------------------
nrf2_tf_tab <- tf_tab %>%
  filter(source_genesymbol %in% c("NFE2L2")) %>%
  select(source_genesymbol, target_genesymbol)

# add some canonical targets from literature ------------------------------
nrf_tf_manual <- rbind(nrf_tf_manual1, nrf_tf_manual2)
nrf_tf_manual <- nrf_tf_manual %>%
  mutate(source_genesymbol = "NFE2L2") %>%
  select(source_genesymbol, target_genesymbol)
## merge
nrf2_tf_tab <- rbind(nrf_tf_manual, nrf2_tf_tab) %>%
  unique()

# write table -------------------------------------------------------------
write.table(x = nrf2_tf_tab, file = paste0(dir_out, "NRF2_Target_Genes.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
