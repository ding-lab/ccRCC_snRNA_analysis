# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input meta data
id_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210322.v1/meta_data.20210322.v1.tsv")
## input ORA enrichment result so as to subset the gene sets later
ora_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/ora_msigdb_tumor_manualsubcluster_up_degs/20210413.v1/ORA_Results.tsv")
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210421.v3/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

# preprocess for subsetting--------------------------------------------------------------
count_geneset_df <- ora_df %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(Description, easy_id) %>%
  unique() %>%
  dplyr::select(Description) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::rename(Description = ".") %>%
  arrange(desc(Freq)) %>%
  mutate(scoregroup_name = paste0(gsub(x = Description, pattern = "HALLMARK_", replacement = ""), "_Score")) %>%
  head(15)
scoregroups_process <- count_geneset_df$scoregroup_name

# make plot data ----------------------------------------------------------
# easyid_tmp <- "C3L-00079-T1"
easyid_tmp <- "C3L-01313-T1"
plotdata_pre_df <- enrich_df %>%
  mutate(cluster_name2 = gsub(x = cluster_name, pattern = "\\.", replacement = "-")) %>%
  filter(grepl(x = cluster_name2, pattern = easyid_tmp))
plotdata_df <- plotdata_pre_df[, scoregroups_process]
rownames(plotdata_df) <- plotdata_pre_df$cluster_name2
plotdata_df <- rbind(rep(100,ncol(plotdata_df)) , rep(-100,ncol(plotdata_df)) , plotdata_df)

# plot --------------------------------------------------------------------
radarchart(plotdata_df)

