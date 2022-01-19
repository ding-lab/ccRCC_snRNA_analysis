# Yige Wu @WashU Sep 2020
## reference: https://www.wikipathways.org/index.php/Pathway:WP4239#nogo2
## reference: https://www.nature.com/articles/s41580-018-0080-4
## Ref: https://pubmed.ncbi.nlm.nih.gov/16567498/

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(clusterProfiler)
library(biomaRt)
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
## input collagens
collagens_df <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/Collagens.csv", skip = 1)
## input keratin
keratins_df <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/Keratins.csv", skip = 1)

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

# manually add EMT genes ------------------------------------------------------------
genes_emt_df <- data.frame(hgnc_symbol = c("VIM", "TWIST1", "ITGB1", "ITGB3", "MMP3", "SOX10", "GCS",
                                           "ITGA3", 
                                           emt_wikipathway_df$hgnc_symbol),
                           gene_function = "Mesenchymal",
                           gene_group = NA,
                           PMID = NA)
genes_emt_df <- rbind(genes_emt_df,
                      data.frame(hgnc_symbol = collagens_df$`Approved symbol`,
                                 gene_function = "Unknown",
                                 gene_group = NA,
                                 PMID = NA))
genes_emt_df <- rbind(genes_emt_df,
                      data.frame(hgnc_symbol = keratins_df$`Approved symbol`,
                                 gene_function = "Epithelial",
                                 gene_group = "Keratin",
                                 PMID = NA))
genes_emt_df <- unique(genes_emt_df)
genes_emt_df$gene_group[genes_emt_df$hgnc_symbol %in% collagens_df$`Approved symbol`] <- "Collagen"
genes_emt_df$PMID[genes_emt_df$hgnc_symbol %in% c("ITGA5")] <- "31142737"
genes_emt_df$PMID[genes_emt_df$hgnc_symbol %in% c("KRT18")] <- "15897741, 15897741"
genes_emt_df$PMID[genes_emt_df$hgnc_symbol %in% c("KRT19")] <- "10792488"
genes_emt_df$PMID[genes_emt_df$hgnc_symbol %in% c("KRT8")] <- "10792488, 15897741"
genes_emt_df$PMID[genes_emt_df$hgnc_symbol %in% c("VIM")] <- "10792488, 15897741"
genes_emt_df$PMID[genes_emt_df$hgnc_symbol %in% c("ITGA3")] <- "10792488, 23786209"
## annnotate the gene function to epithelial
tmp <- genes_emt_df$gene_function
### Epithelial cells are normally held together by lateral cell–cell junctions (tight junctions, adherens junctions, gap junctions and desmosomes). They exhibit apical–basal polarity and interact with the underlying basement membrane via hemidesmosomes and α6β4 integrins.
tmp[genes_emt_df$hgnc_symbol %in% c("CDH1", "ITGA6", "ITGB4", "OCLN", "DSP", "JUP", "PKP1", "PKP2", "CRB3", "MPP5", "TJP1")]  <- "Epithelial"
tmp[grepl(x = genes_emt_df$hgnc_symbol, pattern = "CLDN")] <- "Epithelial"
tmp[genes_emt_df$hgnc_symbol %in% collagens_df$`Approved symbol`]  <- "Unknown"
genes_emt_df$gene_function <- tmp
## annnotate the gene function to epithelial
tmp <- genes_emt_df$gene_group
tmp[grepl(x = genes_emt_df$hgnc_symbol, pattern = "CLDN")] <- "Claudin"
genes_emt_df$gene_group <- tmp
## annotate the key EMT genes
## Ref: https://pubmed.ncbi.nlm.nih.gov/16567498/
genes_emt_df <- genes_emt_df %>%
  mutate(Key_Mesenchymal_Genes = ifelse(hgnc_symbol %in% c("VIM", "CDH2", "FOXC2", "SNAI1", "SNAI2", "TWIST1", "FN1", "ITGA4", "ITGB6", "MMP2", "MMP3", "MMP9", "SOX10", "GSC"), 
                                TRUE, 
                                FALSE)) %>%
mutate(Key_Epithelial_Genes = ifelse(hgnc_symbol %in% c("CDH1", "DSP", "OCLN"), 
                              TRUE, 
                              FALSE))
genes_emt_df <- unique(genes_emt_df)

# write ouputs ------------------------------------------------------------
file2write <- paste0(dir_out, "EMT_Genes.", run_id, ".tsv")
write.table(x = genes_emt_df, file = file2write, quote = F, sep = "\t", row.names = F)



