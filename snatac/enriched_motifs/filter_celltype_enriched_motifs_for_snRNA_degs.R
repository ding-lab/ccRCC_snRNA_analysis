# Yige Wu @WashU Sep 2020

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
## input the united enriched motifs
motifs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/unite_celltype_enriched_motifs/20200908.v1/Enriched_Motifs.chromvar.MergedObj.byCell_type.20200908.v1.tsv")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200904.v1/31Aliquot.Barcode2CellType.20200904.v1.tsv", data.table = F)
## specify top n enriched tf by cell group
n_top <- 100

# get the cell groups for each cell type ----------------------------------
celltype2cellgroup_df <- barcode2celltype_df %>%
  select(Cell_type.detailed, Cell_type.shorter, Cell_group.detailed) %>%
  unique() %>%
  mutate(Cell_type.filename = gsub(x = Cell_type.shorter, pattern = "\\/", replacement = "_"))

# filter motifs by cell group ---------------------------------------------
## map cell group
motifs_df$Cell_group.detailed <- mapvalues(x = motifs_df$Cell_type.filename, from = celltype2cellgroup_df$Cell_type.filename, to = celltype2cellgroup_df$Cell_group.detailed)
### change detailed normal epithelial cell type to cell group
idx_change <- (motifs_df$Cell_group.detailed %in% celltype2cellgroup_df$Cell_type.detailed)
motifs_df$Cell_group.detailed[idx_change] <- mapvalues(x = motifs_df$Cell_type.filename[idx_change], from = celltype2cellgroup_df$Cell_type.detailed, to = celltype2cellgroup_df$Cell_group.detailed)
unique(motifs_df$Cell_group.detailed)
## filter
motifs_filtered_df <- motifs_df %>%
  filter(p_val_adj < 0.01) %>%
  filter(pct_diff > 0.1) %>%
  group_by(Cell_group.detailed) %>%
  top_n(wt = avg_logFC, n = n_top)

# make cell type to TF table ----------------------------------------------
motifs_filtered_uniq_df <- motifs_filtered_df %>%
  select(Cell_group.detailed, motif.name) %>%
  unique()
idx_rep <- sapply(X = motifs_filtered_uniq_df$motif.name, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "(var.2)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  length_genes <- length(vec_genes)
  return(length_genes)
})
celltype2tf_df <- motifs_filtered_uniq_df[rep(1:nrow(motifs_filtered_uniq_df), idx_rep),]
genesymbols_tf <- sapply(X = motifs_filtered_uniq_df$motif.name, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "\\(var.2\\)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  return(vec_genes)
})
celltype2tf_df$tf_genesymbol <- unlist(genesymbols_tf)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Enriched_Motifs.", "Top", n_top, "byCellGroup.", "chromvar.MergedObj.byCell_type.",  run_id, ".tsv")
write.table(x = motifs_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "CellType2TF.", "Top", n_top, "byCellGroup.", "Enriched_Motifs.chromvar.MergedObj.byCell_type.", run_id,".tsv")
write.table(x = celltype2tf_df, file = file2write, quote = F, sep = "\t", row.names = F)

