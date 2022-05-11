# Yige Wu @WashU Dec 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "tximport"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}


# input dependencies ------------------------------------------------------
dir_parent <- "./Resources/Bulk_Processed_Data/mRNA/Kallisto/Outputs/Knockout_Cell_Lines/"
paths_abundance <- list.files(path = dir_parent, recursive = T, full.names = T)
paths_abundance <- paths_abundance[grepl(x = paths_abundance, pattern = "abundance")]

# input files -------------------------------------------------------------
abundance_df <- NULL
sample_ids_vec <- NULL
for (i in 1:length(paths_abundance)) {
  path_abundance <- paths_abundance[i]
  sample_id <- str_split(string = path_abundance, pattern = "\\/")[[1]]
  sample_id <- sample_id[length(sample_id) - 1]
  sample_id <- str_split(string = sample_id, pattern = "\\.")[[1]][1]
  sample_ids_vec <- c(sample_ids_vec, sample_id)
  ab_df <- fread(input = path_abundance)
  if (i == 1) {
    id_df <- str_split_fixed(string = ab_df$target_id, pattern = "\\|", n = 9)
    id_df <- data.frame(id_df)
    colnames(id_df) <- c("ensembl_transcript_id", "ensembl_gene_id", "havana_gene_id", "havana_transcript_id", "transcript_name", "gene_name", "bp", "transcript_type")
    abundance_df <- id_df
    abundance_df[, sample_id] <- ab_df$tpm
  } else {
    abundance_df[, sample_id] <- ab_df$tpm
  }
}

# summarize by gene -------------------------------------------------------
# sample_ids <- str_split_fixed(string = paths_abundance, pattern = "\\/", n = 9)[,8]
# sample_ids <- str_split_fixed(string = sample_ids, pattern = "\\.", n = 2)[,1]
names(paths_abundance) <- sample_ids_vec
tx2gene <- id_df[, c("ensembl_transcript_id", "ensembl_gene_id")]

txi.kallisto.tsv <- tximport(paths_abundance, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
abundance_bygene_df <- data.frame(txi.kallisto.tsv$abundance)
abundance_bygene_df$ensembl_gene_id <- rownames(abundance_bygene_df)
abundance_bygene_df$gene_name <- mapvalues(x = abundance_bygene_df$ensembl_gene_id, from = id_df$ensembl_gene_id, to = as.vector(id_df$gene_name))

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "Knockout_celllines.Kallisto.gene_level.TPM.", run_id, ".tsv")
write.table(x = abundance_bygene_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "Knockout_celllines.Kallisto.transcript_level.TPM.", run_id, ".tsv")
write.table(x = abundance_df, file = file2write, sep = "\t", row.names = F, quote = F)

