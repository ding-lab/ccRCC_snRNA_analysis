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
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210805.v1/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")
emt_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores_wgeneexpression_wPTclusters/20210929.v1/MSigDB.Hallmark.tsv")
## input score pre-calculated
epi_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_epithelial_scores_EMTmoduledown_wPT/20210908.v1/EpithelialScore.tsv")
## add PT segment scores
s12_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_PTS12_scores_wPT/20210908.v1/PTS12Score.tsv")
s3_scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_PTS3_scores_wPT/20210908.v1/PTS12Score.tsv")

# make column annotation --------------------------------------------------
cluster_anno_df <- merge(x = emt_scores_df %>%
                           rename(EMT_score = EPITHELIAL_MESENCHYMAL_TRANSITION_Score) %>%
                           select(cluster_name, EMT_score), 
                         y = epi_scores_df %>%
                           rename(epithelial_score = score), by = c("cluster_name"), all = T)
cutoff_epithelial_strong <- quantile(x = cluster_anno_df$epithelial_score[!(grepl(pattern = "PT", cluster_anno_df$cluster_name))], 0.7)
cutoff_epithelial_intermediate <- quantile(x = cluster_anno_df$epithelial_score[!(grepl(pattern = "PT", cluster_anno_df$cluster_name))], 0.4)
cutoff_epithelial_weak <- quantile(x = cluster_anno_df$epithelial_score[!(grepl(pattern = "PT", cluster_anno_df$cluster_name))], 0.2)
# cluster_anno_df <- cluster_anno_df %>%
#   mutate(epithelial_group = ifelse(grepl(pattern = "PT", cluster_anno_df$cluster_name), "PT", 
#                                    ifelse(epithelial_score >= cutoff_epithelial_strong, "Epithelial-strong",
#                                           ifelse(epithelial_score >= cutoff_epithelial_intermediate, "Epithellal-intermediate",
#                                                  ifelse(cluster_name %in% enrich_df$cluster_name[enrich_df$EMT] & epithelial_score <= cutoff_epithelial_weak, "EMT", "Epithelial-weak")))))
cluster_anno_df$PT_S12_score <- mapvalues(x = cluster_anno_df$cluster_name, from = s12_scores_df$cluster_name, to = as.vector(s12_scores_df$score)); cluster_anno_df$PT_S12_score <- as.numeric(cluster_anno_df$PT_S12_score)
cluster_anno_df$PT_S3_score <- mapvalues(x = cluster_anno_df$cluster_name, from = s3_scores_df$cluster_name, to = as.vector(s3_scores_df$score)); cluster_anno_df$PT_S3_score <- as.numeric(cluster_anno_df$PT_S3_score)
cutoff_pt_s12 <- quantile(cluster_anno_df$PT_S12_score, 0.9)
cutoff_pt_s3 <- quantile(cluster_anno_df$PT_S3_score, 0.9)

cluster_anno_df <- cluster_anno_df %>%
  mutate(epithelial_group = ifelse(grepl(pattern = "PT", cluster_anno_df$cluster_name), "PT", 
                                   ifelse(epithelial_score >= cutoff_epithelial_strong, "Epi-H",
                                          ifelse(epithelial_score >= cutoff_epithelial_intermediate, "Epi-M",
                                                 ifelse(cluster_name %in% enrich_df$cluster_name[enrich_df$EMT] & epithelial_score <= cutoff_epithelial_weak, "EMT", "Epi-L"))))) %>%
  mutate(s123_group = ifelse(PT_S12_score >= cutoff_pt_s12,
                             ifelse(PT_S3_score >=  cutoff_pt_s3, "mixed S1/2/3 identity", "S1/2 enriched"), 
                             ifelse(PT_S3_score >=  cutoff_pt_s3, "S3 enriched", "weak segmental identity"))) %>%
  mutate(cluster_name.figure = gsub(x = cluster_name, pattern = "\\.", replacement = "")) %>%
  mutate(cluster_name.figure = gsub(x = cluster_name.figure, pattern = "C3L0|C3N0", replacement = "P"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Tumorcluster_EpithelialGroup.", run_id, ".tsv")
write.table(x = cluster_anno_df, file = file2write, sep = "\t", row.names = F, quote = F)
