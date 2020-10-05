# Yige Wu @WashU Sep 2020

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
## input summary for cell-cell interactin
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/other/summarize_top_exp_celltype_pergene/20200929.v1/cellphonedb.summary_across_celltypes_by_pair.min5.singlegeneranked.tsv")
## input complex name to gene symbol
complex2gene_list <- readRDS(file = "./Resources/Analysis_Results/cell_cell_interaction/other/summarize_top_exp_celltype_pergene/20200929.v1/cellphonedb.complex2gene.RDS")
## input general druggable genes
druggenes_df <- fread(data.table = F, input = "../ccRCC_Drug/Resources/Knowledge/Gene_Lists/approved_target_ids_all.csv")

# annotate ----------------------------------------------------------------
complex2druggable_list <- sapply(complex2gene_list, function(x, gene_vec) {
  x_hit <- x[x %in% gene_vec]
  if (length(x_hit) <= 1) {
    return(x_hit)
  } else {
    x_hit <- paste0(x_hit, collapse = "_")
    return(x_hit)
  }
}, gene_vec = druggenes_df$`Gene Name`)
complex2druggable_vec <- unlist(complex2druggable_list)

summary_df$gene.source.druggable <- complex2druggable_vec[summary_df$gene.source]
summary_df$gene.target.druggable <- complex2druggable_vec[summary_df$gene.target]

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
write.table(x = summary_df, file = file2write, quote = F, sep = "\t", row.names = F)
