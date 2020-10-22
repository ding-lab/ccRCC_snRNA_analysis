# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(monocle)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
## input paths to the monocle objects
obj_monocle <- readRDS(file = "./Resources/snRNA_Processed_Data/Monocle/outputs/CellTypeVer.20200828.v1/C3L-00079_PooledPT_SelfFib_ByCellType/combined_subset_pseudotime_qval_1e-10.rds")
## input meta data
barcode2manualsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_C3L-00079_tumorclusterid/20201008.v1/TumorCellReclustered.BarcodeInfo.20201008.v1tsv", data.table = F)

# prepare data ---------------------------------------------------------
pdata_df <- as.data.frame(pData(obj_monocle))
pdata_df <- pdata_df %>%
  mutate(id_cell = paste0(orig.ident, "_", original_barcode))
barcode2manualsubcluster_df <- barcode2manualsubcluster_df %>%
  mutate(id_cell = paste0(orig.ident, "_", barcode_tumorcellreclustered)) %>%
  mutate(Name_Cluster = paste0("C", (id_manual_cluster+1)))
## map tumor subcluster id
pdata_df$Name_Cluster <- mapvalues(x = pdata_df$id_cell, from = barcode2manualsubcluster_df$id_cell, to = as.vector(barcode2manualsubcluster_df$Name_Cluster))
pdata_df$Name_Cluster[pdata_df$Name_Cluster == pdata_df$id_cell] <- "Non-tumor"

# count -------------------------------------------------------------------
count_tumorcluster_perstate_df <- pdata_df %>%
  filter(Name_Cluster != "Non-tumor") %>%
  select(Name_Cluster, State) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0)
sum_cells_perstate_df <- pdata_df %>%
  filter(Name_Cluster != "Non-tumor") %>%
  select(State) %>%
  table() %>%
  as.data.frame() %>%
  rename(State = ".")
count_tumorcluster_perstate_df$Total_cells_perstate <- mapvalues(x = count_tumorcluster_perstate_df$State, from = sum_cells_perstate_df$State, to = as.vector(sum_cells_perstate_df$Freq))
count_tumorcluster_perstate_df$Total_cells_perstate <- as.numeric(as.vector(count_tumorcluster_perstate_df$Total_cells_perstate))
count_tumorcluster_perstate_df <- count_tumorcluster_perstate_df %>%
  mutate(Frac_bycluster = (Freq/Total_cells_perstate))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "count_tumorcluster_perstate.", "tsv")
write.table(x = count_tumorcluster_perstate_df, file = file2write, quote = F, sep = "\t", row.names = F)


