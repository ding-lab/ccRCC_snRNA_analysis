# Yige Wu @WashU Feb 2022

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
# dir_out <- paste0(makeOutDir(), run_id, "/")
# dir.create(dir_out)

# input dependencies ------------------------------------------------------
cellnumber_percluster_orig_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/count/count_cellnumber_per_manual_cluster_rm_doublet/20210805.v1/CellNumberPerTumorManualCluster.20210805.v1.tsv")
barcode2manualcluster_ds2000_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id_downsample2000cells/20220228.v1/Barcode2TumorSubclusterId.20220228.v1.tsv")

# preprocess --------------------------------------------------------------
count_tumorcluster_orig_df <- cellnumber_percluster_orig_df %>%
  filter(Freq >= 50) %>% ## 95 tumor clusters with over 50 cells after doublet removal
  select(easy_id) %>%
  table() %>%
  as.data.frame() %>%
  rename(easy_id = ".") %>%
  rename(Freq.orig = Freq)

# compare the original vs. downsample 2000 cell version -------------------
cellnumber_percluster_ds2000_df <- barcode2manualcluster_ds2000_df %>%
  # filter(!is.na(id_manual_cluster_w0)) %>%
  group_by(Cluster_Name) %>%
  summarise(Freq = n())

count_tumorcluster_ds2000_df <- barcode2manualcluster_ds2000_df %>%
  filter(!is.na(id_manual_cluster_w0)) %>%
  filter(Cluster_Name %in% cellnumber_percluster_ds2000_df$Cluster_Name[cellnumber_percluster_ds2000_df$Freq >= 50]) %>%
  select(easy_id, Cluster_Name) %>%
  unique() %>%
  select(easy_id) %>%
  table() %>%
  as.data.frame() %>%
  rename(easy_id = ".") %>%
  rename(Freq.ds2000 = Freq)
count_merged_df <- merge(x = count_tumorcluster_orig_df, y = count_tumorcluster_ds2000_df, by = c("easy_id"), all.x = T)

summary(count_merged_df$Freq.orig)
summary(count_merged_df$Freq.ds2000)

median(x = count_merged_df$Freq.orig, na.rm = T)
median(x = count_merged_df$Freq.ds2000, na.rm = T)
