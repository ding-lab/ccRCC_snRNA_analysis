# Yige Wu @WashU Apr 2020
## running on local
## unite all the calculated pathway scores

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
## input HIF pathway score
hifscore_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/calculate_pathway_score/calculate_hif_pathway_score/20200407.v1/hif_pathway.avg_pathway_scaled_avgexp.20200407.v1.tsv", data.table = F)
## input PBAF complex score
pbafscore_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/calculate_pathway_score/calculate_pbaf_complex_score/20200407.v1/pbaf_complex.avg_pathway_scaled_avgexp.20200407.v1.tsv", data.table = F)
## input MTOR pathway score
mtorscore_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/calculate_pathway_score/calculate_mtor_pathway_score/20200407.v1/mtor_pathwayavg_pathway_scaled_avgexp.20200407.v1.tsv", data.table = F)
## input ORA result pathway score
ora_pathwayscore_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/calculate_pathway_score/calculate_ora_pathway_score/20200407.v1/ora_pathwayavg_pathway_scaled_avgexp.20200407.v1.tsv", data.table = F)
## input ORA result complex score
ora_complexscore_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/calculate_pathway_score/calculate_ora_complex_score/20200407.v1/ora_complex.avg_pathway_scaled_avgexp.20200407.v1.tsv", data.table = F)

# unite tables ------------------------------------------------------------
pathwayscore_df <- rbind(hifscore_df,
                         pbafscore_df,
                         mtorscore_df,
                         ora_pathwayscore_df,
                         ora_complexscore_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "avg_pathway_scaled_avgexp.", run_id, ".tsv")
write.table(x = pathwayscore_df, file = file2write, quote = F, sep = "\t", row.names = F)
