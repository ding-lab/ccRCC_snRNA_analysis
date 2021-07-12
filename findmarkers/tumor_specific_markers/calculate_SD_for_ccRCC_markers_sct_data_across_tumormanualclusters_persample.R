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
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/avgeexp_tumorcells_sct_data_by_manualcluster_rm_doublets_on_katmai/20210413.v1/AverageExpression_ByManualTumorSubcluster.20210413.v1.tsv", data.table = F)
## input cell number per cluster
cellnumber_percluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/count_cellnumber_per_manual_cluster_rm_doublet/20210413.v1/CellNumberPerTumorManualCluster.20210413.v1.tsv")
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210712.v1/ccRCC_markers.Surface.20210712.v1.tsv")

# specify genes to filter -------------------------------------------------
## identify clusters with sufficient cell number
cluster_pass_df <- cellnumber_percluster_df %>%
  filter(Freq >= 50)%>%
  mutate(colname_exp = gsub(x = id_cluster_uniq,pattern = "\\-", replacement = "."))

## add name for the marker groups
genes2filter <- markers_df$Gene

# format expression data --------------------------------------------------
exp_filtered_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycluster_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  dplyr::filter(id_bycluster_byaliquot %in% cluster_pass_df$colname_exp) %>%
  mutate(easyid_column = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(cluster_name = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,2])
## filter out non-tumor and NA tumor cluster
## filter out clusters with <50 cells
exp_filtered_long_df <- exp_filtered_long_df %>%
  filter(!(cluster_name %in% c("", "CNA"))) %>%
  filter(id_bycluster_byaliquot %in% cluster_pass_df$colname_exp)

# filter genes based on variation -----------------------------------------
sd_bygene_bysample_df <- exp_filtered_long_df %>%
  group_by(easyid_column, V1) %>%
  summarise(sd_bysample = sd(x = value))

sd_bygene_mean_df <- sd_bygene_bysample_df %>%
  group_by(V1) %>%
  summarize(sd_bysample_mean = mean(sd_bysample, na.rm = T))

# write outputs -----------------------------------------------------------
file2write <- paste0(dir_out, "SD_for_ccRCC_markers_across_tumormanualclusters_persample", ".tsv")
write.table(x = sd_bygene_bysample_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "Mean_SD_for_ccRCC_markers_across_tumormanualclusters_persample", ".tsv")
write.table(x = sd_bygene_mean_df, file = file2write, sep = "\t", row.names = F, quote = F)

