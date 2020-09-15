# Yige Wu @WashU Sep 2020
## reference: https://www.nature.com/articles/s41585-018-0052-7#Sec27
## reference: https://www.nature.com/articles/nrneph.2016.133#Sec3

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(clusterProfiler)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input HIF targets
hif_tf_df <- fread(input = "./Resources/Analysis_Results/dependencies/write_hif_targets/20200428.v1/HIF_Target_Genes.20200428.v1.tsv", data.table = F)
## input protein protein interaction table
ppi_df <- fread(data.table = F, input = "./Resources/Knowledge/PPI/protein_pair_table_v2.txt")
## input kegg 
gene2pathway_df <- fread(data.table = F, input = "./Resources/Knowledge/Databases/CTD_genes_pathways.tsv.gz", col.names = c("GeneSymbol", "GeneID", "PathwayName", "PathwayID"))

# merge genes related to VHL deficiency----------------------------------------------------
## HIF dependent
genes_hif_related_df <- rbind(hif_tf_df %>%
                                mutate(relation_source2target = "TF"), 
                              ppi_df %>%
                                filter(GENE %in% c("HIF1A", "EPAS1")) %>%
                                filter(SUB_GENE.is_complex_partner) %>%
                                rename(source_genesymbol = GENE) %>%
                                rename(target_genesymbol = SUB_GENE) %>%
                                mutate(target_genefunction = NA) %>%
                                mutate(relation_source2target = "complex partner") %>%
                                select(source_genesymbol, target_genesymbol, target_genefunction, relation_source2target))
## others
genes_vhl_related_df <- unique(ppi_df %>%
                                 filter(GENE %in% c("VHL")) %>%
                                 filter(SUB_GENE.is_complex_partner) %>%
                                 rename(source_genesymbol = GENE) %>%
                                 rename(target_genesymbol = SUB_GENE) %>%
                                 mutate(target_genefunction = NA) %>%
                                 mutate(relation_source2target = "complex partner") %>%
                                 select(source_genesymbol, target_genesymbol, target_genefunction, relation_source2target))

## merge
genes_vhl_hif_df <- rbind(genes_hif_related_df, genes_vhl_related_df)
genes_vhl_hif_df$pathway_name <- "VHL-HIF"

# merge genes related to the chromatin remodeling SMGs --------------------
genes_epigenetic_smgs_related_df <- rbind(ppi_df %>%
                                            filter(GENE %in% c("PBRM1", "SETD2", "BAP1", "KDM5C", "KDM6A")) %>%
                                            filter(SUB_GENE.is_complex_partner) %>%
                                            rename(source_genesymbol = GENE) %>%
                                            rename(target_genesymbol = SUB_GENE) %>%
                                            mutate(target_genefunction = ifelse(source_genesymbol == "PBRM1" & target_genesymbol %in% pbaf_genes, "PBAF complex", NA)) %>%
                                            mutate(relation_source2target = "complex partner") %>%
                                            select(source_genesymbol, target_genesymbol, target_genefunction, relation_source2target),
                                          data.frame(source_genesymbol = c("PBRM1", "SETD2", "BAP1", "KDM5C", "KDM6C"),
                                                     target_genesymbol = c("PBRM1", "SETD2", "BAP1", "KDM5C", "KDM6C"),
                                                     target_genefunction = "SMG",
                                                     relation_source2target = "self"))
genes_epigenetic_smgs_related_df <- rbind(genes_epigenetic_smgs_related_df,
                                          data.frame(source_genesymbol = "BAP1",
                                                     target_genesymbol = c("EZH1", "EZH2"),
                                                     target_genefunction = "Polycomb repressive complex 2",
                                                     relation_source2target = "indirect association"))
genes_epigenetic_smgs_related_df <- rbind(genes_epigenetic_smgs_related_df,
                                          data.frame(source_genesymbol = "SETD2",
                                                     target_genesymbol = c("POLR2A", "RBP1", 
                                                                           "PHF1", "PSIP1", "MORF4L1"),
                                                     target_genefunction = c("transcriptional regulation", "transcriptional regulation", 
                                                                             "DNA DSB repair", "DNA DSB repair", "DNA DSB repair"),
                                                     relation_source2target = c("bind", "bind", "target gene bind H3K36me3", "target gene read H3K36me3", "target gene read H3K36me3")))
genes_epigenetic_smgs_related_df$pathway_name <- "Epigenetic machinary"

# merge genes related to PI3K-AKT-MTOR pathway --------------------
genes_pi3kmtor_df <- data.frame(source_genesymbol = NA,
                                target_genesymbol = gene2pathway_df$GeneSymbol[gene2pathway_df$PathwayName == "mTOR signaling pathway"],
                                target_genefunction = NA,
                                relation_source2target = "mTOR signaling pathway",
                                pathway_name = "PI3K-AKT-mTOR Signaling")
genes_rtk_df <- data.frame(source_genesymbol = NA,
                           target_genesymbol = gene2pathway_df$GeneSymbol[gene2pathway_df$PathwayName == "EGFR tyrosine kinase inhibitor resistance"],
                           target_genefunction = NA,
                           relation_source2target = "EGFR tyrosine kinase inhibitor resistance",
                           pathway_name = "PI3K-AKT-mTOR Signaling")
genes_pi3kmtor_related_df <- rbind(genes_pi3kmtor_df, genes_rtk_df)

# merge genes related to focal adjesion pathway --------------------
genes_adhesion_df <- data.frame(source_genesymbol = NA,
                                target_genesymbol = gene2pathway_df$GeneSymbol[gene2pathway_df$PathwayName == "Focal adhesion"],
                                target_genefunction = NA,
                                relation_source2target = NA,
                                pathway_name = "Focal adhesion")

# merge all pathways ------------------------------------------------------
genes_df <- rbind(genes_vhl_hif_df,
                  genes_epigenetic_smgs_related_df,
                  genes_pi3kmtor_related_df,
                  genes_adhesion_df)

# write table -------------------------------------------------------------
write.table(x = genes_df, file = paste0(dir_out, "ccRCC_Pathogenic_Pathways_Genes.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
