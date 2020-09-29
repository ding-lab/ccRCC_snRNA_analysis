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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cellphonedb output
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/annotate_interactions/annotate_cellphone_out/20200925.v2/cell.phone.res.total.run20200818.filtered.formatted.txt")

# count significant cases, average mean values across samples -----------------------------
summary_df <- cellphone_df %>%
  group_by(pair_cell.types, celltypes.source2target,
           is_ligand_receptor, gene.source, gene.target,  
           paired_cellgroups.general, paired_cellgroups.detailed, Sample_type, 
           Cell_group.source, Cell_group.target, Cell_type.source, Cell_type.target, 
           interacting_pair, variable, 
           is_integrin, secreted) %>%
  dplyr::summarize(number_sig_cases = n(), 
                   avg_sig_mean = mean(value, na.rm = T), 
                   avg_rank_sig_mean.paired_cellgroups.general = mean(rank_sig_mean.paired_cellgroups.general, na.rm = T),
                   avg_rank_sig_mean_same_case = mean(rank_sig_mean_same_case, na.rm = T)) %>%
  arrange(desc(number_sig_cases)) %>%
  rename(paired_celltypes = variable)
  
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cellphonedb.summary_across_samples_by_pair.", run_id, ".tsv")
write.table(x = summary_df, file = file2write, quote = F, sep = "\t", row.names = F)