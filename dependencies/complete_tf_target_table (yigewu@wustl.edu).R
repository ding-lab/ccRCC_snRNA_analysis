# Yige Wu @WashU Oct 2020
## adding known TF-targets to the omnipath TF-target relations

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
## input TF integrations
tf2target_omnipath_df <- fread(input = "./Resources/Knowledge/PPI/Transcriptional/omnipathdb.transcriptional.20200908.txt", data.table = F)
## input known HIF targets
hif2target_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_hif_targets/20201027.v1/HIF_Target_Genes.20201027.v1.tsv")

# merge -------------------------------------------------------------------
tf2target_merged_df <- merge(x = tf2target_omnipath_df, 
                             y = hif2target_df %>%
                               mutate(is_directed = 1), by = c("source_genesymbol", "target_genesymbol", "is_directed"), all = T)

# write output ------------------------------------------------------------
write.table(x = tf2target_merged_df, file = paste0(dir_out, "omnipathdb.transcriptional.plus_manual.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)


