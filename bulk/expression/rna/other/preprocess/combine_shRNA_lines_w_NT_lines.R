# Yige Wu @WashU May 2022
## reference: https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf
## how they deal with count of 0 to log2CPM: https://support.bioconductor.org/p/107719/

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
  "data.table"
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
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input ------------------------------------------------------
counts1_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/Knockdown_Cell_Lines/all.gene_counts.tsv")
counts2_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/Cell_Lines/all.gene_counts.tsv")
counts3_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/BAP1_Cell_Lines/all.gene_counts.tsv")

# preprocess --------------------------------------------------------------
counts2_df <- counts2_df[, 1:15]
colnames_keep <- colnames(counts3_df)
colnames_keep <- c(colnames_keep[1:7], colnames_keep[grepl(pattern = "786|skrc", x = colnames_keep)])
counts3_df <- counts3_df[, colnames_keep]

##
which(is.na(counts1_df$ensembl_gene_id))
which(is.na(counts1_df$entrezgene))
which(is.na(counts1_df$external_gene_name))
which(is.na(counts1_df$gene_biotype))
which(is.na(counts1_df$external_gene_source))
which(is.na(counts1_df$transcript_count))
which(is.na(counts1_df$description))
## only entrezgene has NAs
counts1_df$entrezgene <- as.character(counts1_df$entrezgene); counts1_df$entrezgene[is.na(counts1_df$entrezgene)] <- ""
counts2_df$entrezgene <- as.character(counts2_df$entrezgene); counts2_df$entrezgene[is.na(counts2_df$entrezgene)] <- ""
counts3_df$entrezgene <- as.character(counts3_df$entrezgene); counts3_df$entrezgene[is.na(counts3_df$entrezgene)] <- ""

## merge
counts_merged_df <- merge(x = counts1_df %>%
                            dplyr::select(-external_gene_source), 
                          y = counts2_df %>%
                            dplyr::select(-external_gene_source), by = c("ensembl_gene_id", "entrezgene", "external_gene_name"), all = T, suffixes = c(".1", ".2"))
counts_merged_df <- merge(x = counts_merged_df, 
                          y = counts3_df %>%
                            dplyr::select(-external_gene_source), by = c("ensembl_gene_id", "entrezgene", "external_gene_name"), all = T)

which(counts_merged_df$description.x !=  counts_merged_df$description.y)
colnames_merged <- colnames(counts_merged_df)
counts_merged_df <- counts_merged_df[, c(colnames_merged[!grepl(pattern = "sample", x = colnames_merged)], colnames_merged[grepl(pattern = "sample", x = colnames_merged)])]


## write output
file2write <- paste0(dir_out, "Cell_Lines.gene_counts.", run_id, ".tsv")
write.table(x = counts_merged_df, file = file2write, row.names = F, quote = F, sep = "\t")

