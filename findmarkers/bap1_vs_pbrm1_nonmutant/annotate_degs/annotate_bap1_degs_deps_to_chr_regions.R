# Yige Wu @WashU Jun 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(biomaRt)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/unite_BAP1_snRNA_bulkRNA_protein_DEGs/20210616.v1/BAP1_DEGs.United.snRNA.bulkRNA.Protein.20210616.v1.tsv")
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
# filters = listFilters(ensembl)
# attributes_ens <- listAttributes(ensembl)

# retrieve genomic region -------------------------------------------------
genes2convert <- deg_df$genesymbol_deg
genesymbol2region_df <- getBM(attributes=c('hgnc_symbol', 'chromosome_name', "start_position", "end_position", "strand", "band"), 
                                filters = 'hgnc_symbol', 
                                values = genes2convert, 
                                mart = ensembl)
genesymbol2region_df <- genesymbol2region_df %>%
  arrange(hgnc_symbol, desc(band))
genesymbol2region_filtered_df <- genesymbol2region_df[!duplicated(genesymbol2region_df$hgnc_symbol),]
deg_annotated_df <- merge(x = deg_df, y = genesymbol2region_filtered_df, by.x = c("genesymbol_deg"), by.y = c("hgnc_symbol"), all.x = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_vs_PBRM1_NonMutants.DEGs.Chromosome_Regions.", run_id, ".tsv")
write.table(x = genesymbol2region_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_vs_PBRM1_NonMutants.DEGs.Chromosome_Regions_Annotated.", run_id, ".tsv")
write.table(x = deg_annotated_df, file = file2write, quote = F, sep = "\t", row.names = F)
