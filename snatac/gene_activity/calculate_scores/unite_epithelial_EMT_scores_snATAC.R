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
# ## input by cluster enrichment assignment
emt_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/gene_activity/calculate_scores/calculate_msigdb_geneset_scores_wgeneactivity/20210921.v1/MSigDB.Hallmark.tsv")
epi_group_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_epithelial_group_bytumorcluster/20210929.v1/Tumorcluster_EpithelialGroup.20210929.v1.tsv")
## input score pre-calculated
epi_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/gene_activity/calculate_scores/calculate_Epithelial_scores_EMTmoduledown_wPT_wgeneactivity/20210921.v1/EpithelialScore.tsv")
s12_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/gene_activity/calculate_scores/calculate_PTS12_scores_wPT_wgeneactivity/20210921.v1/PTS12Score.tsv")
s3_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/gene_activity/calculate_scores/calculate_PTS3_scores_wPT_wgeneactivity/20210921.v1/PTS3Score.tsv")

# unite -------------------------------------------------------------------
cluster_anno_df <- epi_scores_df
cluster_anno_df$PT_S12_score <- mapvalues(x = cluster_anno_df$cluster_name, from = s12_scores_df$cluster_name, to = as.vector(s12_scores_df$score)); cluster_anno_df$PT_S12_score <- as.numeric(cluster_anno_df$PT_S12_score)
cluster_anno_df$PT_S3_score <- mapvalues(x = cluster_anno_df$cluster_name, from = s3_scores_df$cluster_name, to = as.vector(s3_scores_df$score)); cluster_anno_df$PT_S3_score <- as.numeric(cluster_anno_df$PT_S3_score)
cluster_anno_df$EMT_score <- mapvalues(x = cluster_anno_df$cluster_name, from = emt_scores_df$cluster_name, to = as.vector(emt_scores_df$EPITHELIAL_MESENCHYMAL_TRANSITION_Score)); cluster_anno_df$EMT_score <- as.numeric(cluster_anno_df$EMT_score)
cluster_anno_df$epithelial_group <- mapvalues(x = cluster_anno_df$cluster_name, from = epi_group_df$cluster_name, to = as.vector(epi_group_df$epithelial_group));
## assign S12/S3 signature
cutoff_pt_s12 <- quantile(cluster_anno_df$PT_S12_score, 0.9)
cutoff_pt_s3 <- quantile(cluster_anno_df$PT_S3_score, 0.9)
cluster_anno_df <- cluster_anno_df %>%
  mutate(s123_group = ifelse(PT_S12_score >= cutoff_pt_s12,
                             ifelse(PT_S3_score >=  cutoff_pt_s3, "mixed S1/2/3 identity", "S1/2 enriched"), 
                             ifelse(PT_S3_score >=  cutoff_pt_s3, "S3 enriched", "weak segmental identity"))) %>%
  mutate(cluster_name.figure = gsub(x = cluster_name, pattern = "\\.", replacement = "")) %>%
  mutate(cluster_name.figure = gsub(x = cluster_name.figure, pattern = "C3L0|C3N0", replacement = "P"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Epithelial_EMT_scores_snATACbased_per_tumorcluster.", run_id, ".tsv")
write.table(x = cluster_anno_df, file = file2write, quote = F, sep = "\t", row.names = F)
