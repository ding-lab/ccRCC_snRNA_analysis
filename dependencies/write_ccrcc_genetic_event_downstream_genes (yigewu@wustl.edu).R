# Yige Wu @WashU Feb 2020
## for writing a table for the potential downstream genes for ccRCC genetic alterations

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input HIF targets
hif_tf_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_hif_targets/", data.table = F)
## input MYC targets
myc_tf_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_myc_targets/20200227.v1/MYC_Target_Genes.20200227.v1.tsv", data.table = F)
## input NRF2 targets
nrf2_tf_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_nrf2_targets/20200302.v1/NRF2_Target_Genes.20200302.v1.tsv", data.table = F)
## input TP53 targets
tp53_tf_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/PPI/TF_interactions_TP53_manual.txt", data.table = F)
## input SQSTM1 (5q gain) affected genes
sqstm1_df <- data.frame(source_genesymbol = "SQSTM1",
                        target_genesymbol = c("KEAP1", "NFE2L2", 
                                              "TRAF6", "NFKB1", "MTOR"))
## input MTOR affected genes
mtor_df <- data.frame(source_genesymbol = "MTOR",
                      target_genesymbol = c("SREBP1", "SREBP2"))
### input SREBP target genes
srebp_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_srebp_targets/20200302.v1/SREBP_Target_Genes.20200302.v1.tsv", data.table = F)
# merge ------------------------------------
genetic_alt_downstream_genes <- rbind(tp53_tf_tab[, c("source_genesymbol", "target_genesymbol")],
                                      nrf2_tf_tab, myc_tf_tab, hif_tf_tab,
                                      sqstm1_df, mtor_df, srebp_df)

# write table -------------------------------------------------------------
write.table(x = genetic_alt_downstream_genes, file = paste0(dir_out, "ccRCC_Genetic_Event_Downstream_Genes.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
