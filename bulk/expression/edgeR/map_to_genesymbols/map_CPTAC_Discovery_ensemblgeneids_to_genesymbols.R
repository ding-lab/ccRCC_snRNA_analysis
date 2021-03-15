# Yige Wu @WashU Mar 2021
## some genes only have clone-based ensembl gene id: http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000167046;r=20:62640719-62643304;t=ENST00000370520
## some genes are only in GRCh37.p13 build but not GRCh38: https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000000005;r=X:99839799-99854882

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")
library(biomaRt)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input genes
dgelist_obj <- readRDS(file = './Resources/Analysis_Results/bulk/expression/edgeR/create_DGEList_object/create_DEGList_CPTAC_ccRCC_Discovery_cases/20210312.v1/CPTAC_Discovery_ccRCC_Cases.DEGList.RDS')

# preprocess ------------------------------------------------------------------
genes2convert_df <- data.frame(ensembl_gene_id_version = rownames(dgelist_obj))
rm(dgelist_obj)
genes2convert_df <- genes2convert_df %>%
  mutate(ensembl_gene_id = str_split_fixed(string = ensembl_gene_id_version, pattern = "\\.", n = 2)[,1])
ensembol_gene_ids2convert <- unique(genes2convert_df$ensembl_gene_id)

# run biomart -------------------------------------------------------------
## get biomart
### reference: https://grch37.ensembl.org/info/data/biomart/biomart_r_package.html
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
## retrieve entrezgene_id
mapping_df <- getBM(attributes=c("ensembl_gene_id", 'ensembl_gene_id_version', 'entrezgene_id', 'hgnc_symbol', 'clone_based_ensembl_gene'), 
                                filters = 'ensembl_gene_id', 
                                values = ensembol_gene_ids2convert, 
                                mart = ensembl)
nrow(mapping_df) ## 22391
mapping_nodup_df <- mapping_df[!duplicated(mapping_df$ensembl_gene_id),]

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ensembl_gene_id.mapping.", run_id, ".tsv")
write.table(x = mapping_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ensembl_gene_id.mapping.nodup.", run_id, ".tsv")
write.table(x = mapping_nodup_df, file = file2write, quote = F, sep = "\t", row.names = F)
