# Yige Wu @WashU Sep 2020
## reference of the cellphonedb output: https://www.cellphonedb.org/documentation

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
## input cellphonedb output
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/annotate_interactions/annotate_cellphone_out/20201012.v1/cell.phone.res.total.run20200818.filtered.formatted.txt")
cellphone_sum_by_paircelltypes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/summarize_across_samples_by_pair_and_celltypes/20201012.v1/cellphonedb.summary_across_samples_by_pair.20201012.v1.tsv")


# filter by 79 ------------------------------------------------------------
cellphone_filtered_df1 <- cellphone_df %>%
  filter(gene.source %in% c("WNT5A", "WNT5B")) %>%
  filter(Cell_type.source == "Tumor-like epithelial cells") %>%
  filter(Easy_id %in% c("C3L-00079-T1")) %>%
  arrange(desc(value))
cellphone_filtered_df2 <- cellphone_df %>%
  filter(gene.source %in% c("WNT5A", "WNT5B")) %>%
  filter(Cell_type.source == "Tumor-like epithelial cells") %>%
  filter(Cell_type.target %in% "Tumor-like epithelial cells") %>%
  filter(Easy_id %in% c("C3L-00079-T1")) %>%
  arrange(desc(value))
# There are 3 significant WNT5A interactions ( WNT5A_ROR1, FZD6_WNT5A, WNT5A_ANTXR1) among the tumor-like epithelial cells/transitional cells themselves
cellphone_filtered_df0 <- cellphone_df %>%
  filter(interacting_pair %in% c("FZD6_WNT5A")) %>%
  filter(Easy_id %in% c("C3L-00079-T1")) %>%
  arrange(desc(value))
cellphone_filtered_df01 <- cellphone_df %>%
  filter(interacting_pair %in% c("FZD6_WNT5A")) %>%
  select(celltypes.source2target) %>%
  unique()
dim(cellphone_filtered_df01)

# filter 1200 -------------------------------------------------------------
cellphone_filtered_df1 <- cellphone_df %>%
  filter(gene.source %in% c("WNT5B", "WNT5A") | gene.target %in% c("LRP6")) %>%
  # filter(Cell_type.source == "Tumor cells") %>%
  filter(Easy_id %in% c("C3N-01200-T3")) %>%
  arrange(desc(value))

# filter by the 3 WNT5A interactions --------------------------------------
cellphone_sum_rankby_avg_sig_mean_df <- cellphone_sum_by_paircelltypes_df %>%
  filter(interacting_pair %in% cellphone_filtered_df$interacting_pair) %>%
  arrange(interacting_pair, desc(avg_sig_mean))

cellphone_sum_rankby_avg_sig_mean_filterby_num_sig_case_df <- cellphone_sum_by_paircelltypes_df %>%
  filter(interacting_pair %in% cellphone_filtered_df$interacting_pair) %>%
  filter(number_sig_cases >= 7) %>%
  arrange(interacting_pair, desc(avg_sig_mean))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cellphonedb.summary_across_samples_by_pair.", run_id, ".tsv")
write.table(x = summary_df, file = file2write, quote = F, sep = "\t", row.names = F)
