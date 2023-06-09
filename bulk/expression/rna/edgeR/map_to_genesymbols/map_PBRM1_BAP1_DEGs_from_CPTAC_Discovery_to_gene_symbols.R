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
## input DEGs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/run_deg_analysis/run_PBRM1_BAP1_deg_on_cptac_ccRCC_discovery_cases/20210312.v1/glmQLFTest.outputtables.tsv")

# examine data ------------------------------------------------------------
length(unique(deg_df$gene_ensembl_id)) # 22589
deg_filtered_df <- deg_df %>%
  filter(PValue < 0.001)
length(unique(deg_filtered_df$gene_ensembl_id)) ## 19464
table(deg_filtered_df$comparison)
# BAP1 mutated vs NAT  BAP1 mutated vs Non-mutants           Non-mutants vs NAT         PBRM1 mutated vs NAT PBRM1 mutated vs Non-mutants 
# 14812                         5053                        15815                        16254                         2659

# preprocess ------------------------------------------------------------------
deg_df <- deg_df %>%
  rename(ensembl_gene_id_version = gene_ensembl_id) %>%
  mutate(ensembl_gene_id = str_split_fixed(string = ensembl_gene_id_version, pattern = "\\.", n = 2)[,1])
ensembol_gene_ids2convert <- unique(deg_df$ensembl_gene_id)

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

# merge -------------------------------------------------------------------
mapping_nodup_df <- mapping_df[!duplicated(mapping_df$ensembl_gene_id),]
deg_df <- merge(x = deg_df, y = mapping_nodup_df, by.x = c("ensembl_gene_id"), by.y = c("ensembl_gene_id"), all.x = T)
deg_df <- deg_df %>%
  filter(!(ensembl_gene_id_version.x %in% c("__alignment_not_unique", "__ambiguous", "__no_feature")))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_BAP1_deg_on_cptac_ccRCC_discovery_cases.w_genesymbols.", run_id, ".tsv")
write.table(x = deg_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ensembl_gene_id.mapping.", run_id, ".tsv")
write.table(x = mapping_df, file = file2write, quote = F, sep = "\t", row.names = F)

