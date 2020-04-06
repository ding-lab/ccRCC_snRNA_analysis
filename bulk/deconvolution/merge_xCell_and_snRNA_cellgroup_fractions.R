# Yige Wu @WashU March 2020
## merge xCell values with snRNA cell group fractions

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input xCell values
xcell_value_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/bulk/deconvolution/summarize_xCell_cell_groups/20200320.v1/xcell_parent_cells_fraction.long.20200320.v1.tsv", data.table = F)
## input snRNA cell group fractions
snRNA_cellgroup_stat_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/other/summarize_snRNA_cell_groups_for_xCell_comparison/20200320.v1/snRNA_cellgroup_stats_for_xCell_comparison.20200320.v1.tsv", data.table = )

# merge -------------------------------------------------------------------
cellgroup_stat_df <- merge(xcell_value_df, snRNA_cellgroup_stat_df, by = c("Aliquot.snRNA", "Cellgroup.xCell"), all = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "merged_xCell_value_and_snRNA_cellgroup_fractions.", run_id, ".tsv")
write.table(x = cellgroup_stat_df, file = file2write, quote = F, sep = "\t", row.names = F)
