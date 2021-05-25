# Yige Wu @WashU Oct 2020

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
## input the DEG-TF matrix
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/overlap_tumor_vs_normal_snRNA_bulkRNA_protein_DEGs/20210511.v1/Tumor_vs_PT_DEGs.Overlap.snRNA.bulkRNA.Protein.20210511.v1.tsv")
## specify enriched motifs
motifs_plot <- c("NFKB2", "NFKB1", "HIF1A", "ARNT::HIF1A", "RBPJ", "MXI1", "ZNF75D", "HSF2", "NEUROD1", "SREBF2", "NEUROG2(var.2)",
                 "KLF15", "NRF1", "SP9", "ZBTB14", "EGR1", "SP3", "TCFL5", "ZNF148", "KLF14", "SP1")

# melt the DEG-TF table ---------------------------------------------------
### get the gene symbols of the TFs
genesymbols_tf <- sapply(X = motifs_plot, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "\\(var.2\\)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  return(vec_genes)
})
deg_filtered_df <- deg_df %>%
  filter(genesymbol_deg %in% unlist(genesymbols_tf))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "EnrichedTFs", ".DifferentialExpression", ".tsv")
write.table(x = deg_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
