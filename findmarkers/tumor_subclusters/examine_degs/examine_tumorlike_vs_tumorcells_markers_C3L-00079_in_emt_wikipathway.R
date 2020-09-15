# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(clusterProfiler)
library(biomaRt)
# library(org.Hs.eg.db)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
## input wikipathway 
wp2gene <- read.gmt("./Resources/Knowledge/Databases/Wikipathways/wikipathways-20200510-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_C3L-00079_tumorlike_vs_tumorcells/20200909.v1/findmarkers_wilcox_tumorlikecells_vs_tumorcells.logfcthreshold0.693147180559945.minpct0.1.mindiffpct0.1.tsv")

# filter wikipathway to the EMT pathway, convert entrez ids to gene symbol ---------------------------------------
emt_wikipathway_df <- wp2gene %>%
  dplyr::filter(name == "Epithelial to mesenchymal transition in colorectal cancer")
genes2convert <- unique(emt_wikipathway_df$gene)
## retrieve entrezgene_id
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'entrezgene_id', 
                                values = genes2convert, 
                                mart = ensembl)
## add entrez ids to the deg table
emt_wikipathway_df$hgnc_symbol <- mapvalues(x = emt_wikipathway_df$gene, from = genesymbol2entrezid_df$entrezgene_id, to = as.vector(genesymbol2entrezid_df$hgnc_symbol))

# filter DEGs -------------------------------------------------------------
deg_emt_df <- deg_df %>%
  dplyr::filter(row_name %in% emt_wikipathway_df$hgnc_symbol)



