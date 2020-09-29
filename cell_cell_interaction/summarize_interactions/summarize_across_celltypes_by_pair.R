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
cellphone_sum_by_paircelltypes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/summarize_across_samples_by_pair_and_celltypes/20200924.v1/cellphonedb.summary_across_samples_by_pair.20200924.v1.tsv")

# specify the least number of significant samples -------------------------
min_number_sig_cases <- 5

# count significant cases, average mean values across samples -----------------------------
nrow(cellphone_sum_by_paircelltypes_df)
cellphone_sum_by_paircelltypes_df <- unique(cellphone_sum_by_paircelltypes_df)
nrow(cellphone_sum_by_paircelltypes_df)
cellphone_sum_by_pair_df <- cellphone_sum_by_paircelltypes_df %>%
  filter(number_sig_cases >= min_number_sig_cases) %>%
  filter(Sample_type == "Tumor") %>%
  dplyr::group_by(interacting_pair, 
                  is_ligand_receptor, gene.source, gene.target,  
                  is_integrin, secreted) %>%
  dplyr::summarize(paired_cellgroups.general = paired_cellgroups.general[which.min(avg_rank_sig_mean_same_case)], 
                   number_sig_cases = number_sig_cases[which.min(avg_rank_sig_mean_same_case)],
                   avg_rank_sig_mean.same_case = avg_rank_sig_mean_same_case[which.min(avg_rank_sig_mean_same_case)],
                   avg_rank_sig_mean.paired_cellgroups.general = avg_rank_sig_mean.paired_cellgroups.general[which.min(avg_rank_sig_mean.paired_cellgroups.general)],
                   avg_sig_mean = avg_sig_mean[which.min(avg_rank_sig_mean_same_case)],
                   paired_celltypes = paired_celltypes[which.min(avg_rank_sig_mean_same_case)], 
                   paired_cellgroups.detailed = paired_cellgroups.detailed[which.min(avg_rank_sig_mean_same_case)], 
                   Cell_group.source = Cell_group.source[which.min(avg_rank_sig_mean_same_case)],
                   Cell_group.target = Cell_group.target[which.min(avg_rank_sig_mean_same_case)],
                   Cell_type.source = Cell_type.source[which.min(avg_rank_sig_mean_same_case)],
                   Cell_type.target = Cell_type.target[which.min(avg_rank_sig_mean_same_case)],
                   interaction_celltypes = pair_cell.types[which.min(avg_rank_sig_mean_same_case)])
# cellphone_sum_by_pair_df <- cellphone_sum_by_paircelltypes_df %>%
#   filter(number_sig_cases >= 5) %>%
#   group_by(interacting_pair) %>%
#   dplyr::summarize(top_paired_celltypes = paired_celltypes[which.max(avg_sig_mean)])

table(cellphone_sum_by_pair_df$paired_cellgroups.general)
table(cellphone_sum_by_pair_df$is_integrin)
nrow(cellphone_sum_by_pair_df)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cellphonedb.summary_across_celltypes_by_pair.min", min_number_sig_cases, ".tsv")
write.table(x = cellphone_sum_by_pair_df, file = file2write, quote = F, sep = "\t", row.names = F)
