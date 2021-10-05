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
cluster_anno_df <- cluster_anno_df %>%
  mutate(epithelial_group = ifelse(grepl(pattern = "PT", cluster_anno_df$cluster_name), "PT", 
                                   ifelse(epithelial_score >= cutoff_epithelial_strong, "Epi-H",
                                          ifelse(epithelial_score >= cutoff_epithelial_intermediate, "Epi-M",
                                                 ifelse(cluster_name %in% enrich_df$cluster_name[enrich_df$EMT] & epithelial_score <= cutoff_epithelial_weak, "EMT", "Epi-L")))))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Tumorcluster_EpithelialGroup.", run_id, ".tsv")
write.table(x = cluster_anno_df, file = file2write, sep = "\t", row.names = F, quote = F)
