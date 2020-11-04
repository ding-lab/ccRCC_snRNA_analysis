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
deg2motif_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_degs_associated_with_tf_with_known_relations/20201022.v1/DEGs_with_TFs_inDARs.Long.tsv")

# filter to rcc pathway genes ---------------------------------------------
deg2motif_filtered_df <- deg2motif_df %>%
  filter(genesymbol_deg %in% c("ENPP3", "EGFR", "PLOD2", "PFKP", "PKM")) %>%
  filter(value == 1) %>%
  arrange(genesymbol_deg)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_Pathway_DEG_Promoter_TFMotifs", ".tsv")
write.table(x = deg2motif_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

