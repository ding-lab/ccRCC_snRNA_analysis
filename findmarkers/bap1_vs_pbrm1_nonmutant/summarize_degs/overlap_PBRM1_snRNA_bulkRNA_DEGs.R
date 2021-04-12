# Yige Wu @WashU Mar 2021

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
deg_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/summarize_PBRM1_vs_BAP1_NonMutant_DEGs/20210408.v1/PBRM1_DEGs.Consistent20210408.v1.tsv")
deg_bulkRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/run_deg_analysis/run_PBRM1_vs_BAP1_NonMutant_deg_on_cptac_ccRCC_discovery_cases/20210408.v1/PBRM1_Mutated_vs_BAP1_NonMutants.glmQLFTest.OutputTables.tsv")

# overlap -----------------------------------------------------------------
## preprocess the bulk DEGs
deg_bulkRNA_sig_df <- deg_bulkRNA_df %>%
  filter(FDR < 0.05) %>%
  mutate(direction.bulkRNA = ifelse(logFC > 0, "Up", "Down")) %>%
  mutate(uniq_id = paste0(hgnc_symbol, "_", direction.bulkRNA))

deg_overlap_df <- merge(x = deg_snRNA_df %>%
                         mutate(uniq_id = paste0(genesymbol_deg, "_", PBRM1_vs_OtherTumor_snRNA)), 
                       y = deg_bulkRNA_sig_df, by = c("uniq_id"))
nrow(deg_overlap_df)
table(deg_overlap_df$direction.bulkRNA)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_DEGs.Overlap.", "snRNA_and_bulkRNA.", run_id, ".tsv")
write.table(x = deg_overlap_df, file = file2write, quote = F, sep = "\t", row.names = F)
