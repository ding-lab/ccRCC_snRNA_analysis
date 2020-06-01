# Yige Wu @WashU May 2020
## run DEG analysis for cells with certain CNV vs CN neutral

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(openxlsx)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_with_cna/findallmarkers_genelevel_expected_cnv_vs_neutral_in_tumorcells/20200518.v1/ExpectedCNV_vs_Neutral..FindAllMarkers.Wilcox..20200518.v1.tsv")

# summarize de gene occurence and unique de genes -------------------------
countaliquot_bygene_df <- markers_df %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cna_gene_symbol, de_gene_symbol) %>%
  summarise(num_aliquot_de = n(), mean_avg_logFC = mean(avg_logFC))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "markers_associated_with_cna.", "p_val_adj.05.", run_id, ".tsv")
write.table(x = countaliquot_bygene_df, file = file2write, quote = F, row.names = F, sep = "\t")

# write output by cna gene for de genes up-regulated in cna samples-----------------------------
list_degs <- list()
for (cna_gene_tmp in unique(countaliquot_bygene_df$cna_gene_symbol)) {
  df2write <- countaliquot_bygene_df %>%
    filter(cna_gene_symbol == cna_gene_tmp) %>%
    filter(mean_avg_logFC > 0) %>%
    arrange(desc(num_aliquot_de))
  list_degs[[cna_gene_tmp]] <- df2write
}
file2write <- paste0(dir_out, "up." , "markers_associated_with_cna.", "p_val_adj.05.", run_id, ".xlsx")
write.xlsx(list_degs, file = file2write)

# write output by cna gene for de genes up-regulated in cna samples-----------------------------
list_degs <- list()
for (cna_gene_tmp in unique(countaliquot_bygene_df$cna_gene_symbol)) {
  df2write <- countaliquot_bygene_df %>%
    filter(cna_gene_symbol == cna_gene_tmp) %>%
    filter(mean_avg_logFC > 0) %>%
    arrange(desc(num_aliquot_de))
  list_degs[[cna_gene_tmp]] <- df2write
}
file2write <- paste0(dir_out, "up." , "markers_associated_with_cna.", "p_val_adj.05.", run_id, ".xlsx")
write.xlsx(list_degs, file = file2write)

# write output by cna gene for de genes down-regulated in cna samples-----------------------------
list_degs <- list()
for (cna_gene_tmp in unique(countaliquot_bygene_df$cna_gene_symbol)) {
  df2write <- countaliquot_bygene_df %>%
    filter(cna_gene_symbol == cna_gene_tmp) %>%
    filter(mean_avg_logFC < 0) %>%
    arrange(desc(num_aliquot_de))
  list_degs[[cna_gene_tmp]] <- df2write
}
file2write <- paste0(dir_out, "down." , "markers_associated_with_cna.", "p_val_adj.05.", run_id, ".xlsx")
write.xlsx(list_degs, file = file2write)

