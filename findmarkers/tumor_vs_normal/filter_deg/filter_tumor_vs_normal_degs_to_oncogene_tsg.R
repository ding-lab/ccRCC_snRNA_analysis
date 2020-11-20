# Yige Wu @WashU Nov 2020

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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_tumor_vs_pt_on_katmai/20200903.v1/findallmarkers_wilcox_tumorcells_vs_pt.20200903.v1.tsv")
## input oncogene and tsgs
cancerdrivers_df <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/Consensus.Genelist.full.txt", header = T)
nrow(cancerdrivers_df)

# filter to rcc pathway genes ---------------------------------------------
deg_filtered_df <- deg_df %>%
  rename(genesymbol_deg = row_name) %>%
  filter(genesymbol_deg %in% cancerdrivers_df$Gene) %>%
  arrange(desc(avg_logFC))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "RCC_Pathway_DEG_Promoter_TFMotifs", ".tsv")
write.table(x = deg2motif_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

