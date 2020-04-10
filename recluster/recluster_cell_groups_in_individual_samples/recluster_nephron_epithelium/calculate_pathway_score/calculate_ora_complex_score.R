# Yige Wu @WashU Apr 2020
## running on local
## for calculating pathway activity score using scaled average expression
## complexes

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
## input scaled average expression
scaled_avgexp_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/other/scale_averageexpression/20200407.v1/averageexpression_tumor_cells_by_manual_subcluster.scaled.20200407.v1.tsv", data.table = F)
## input ORA pathway results
ora_complex_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/findvariablefeatures/findvariablefeatures_tumor_cells/20200310.v1/ORA_results_Complex.tsv", data.table = F)
## input the variable gene list
variable_genes_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/findvariablefeatures/findvariablefeatures_tumor_cells/20200310.v1/findvariablefeatures_tumor_cells.20200310.v1.tsv", data.table = F)

# filter out rows with all NAs in scaled data--------------------------------------------
scaled_avgexp_df <- scaled_avgexp_df[rowSums(!is.na(scaled_avgexp_df)) > 1,]


# process by pathway ------------------------------------------------------
avg_pathway_scaled_avgexp_df <- NULL
pathway_members <- list()
for (path_tmp in ora_complex_df$complex_name) {
  # get the list of eligible HIF pathway members ----------------------------
  pathway_genes_string <- ora_complex_df$members_input_overlap[ora_complex_df$complex_name == path_tmp]
  pathway_genes <- str_split(string = pathway_genes_string, pattern = "; ")[[1]]
  pathway_genes
  
  # average by column/manual subcluster -------------------------------------------------------
  ## filter row by the pathway members
  pathway_scaled_avgexp_df <- scaled_avgexp_df %>%
    filter(V1 %in% pathway_genes)
  ## get matrix-like scaled data
  pathway_scaled_avgexp_mat <- pathway_scaled_avgexp_df %>%
    select(-V1)
  ## average by column
  avg_pathway_scaled_avgexp_mat <- colMeans(pathway_scaled_avgexp_mat)
  ## make long data frame
  avg_pathway_scaled_avgexp_tmp <- data.frame(exp_manual_subcluster_name = str_split_fixed(string = names(avg_pathway_scaled_avgexp_mat), pattern = "\\.", n = 2)[,2],
                                             pathway_score = avg_pathway_scaled_avgexp_mat,
                                             pathway_name = path_tmp)
  avg_pathway_scaled_avgexp_df <- rbind(avg_pathway_scaled_avgexp_tmp, avg_pathway_scaled_avgexp_df)
  
  ## store pathway members
  pathway_members[[path_tmp]] <- pathway_genes
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ora_complex.", "avg_pathway_scaled_avgexp.", run_id, ".tsv")
write.table(x = avg_pathway_scaled_avgexp_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ora_complex.", "pathway_members.", run_id, ".RDS")
saveRDS(object = pathway_members, file = file2write)

