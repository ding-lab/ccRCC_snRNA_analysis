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
## input HIF targets
hif_tf_df <- fread(input = "./Resources/Analysis_Results/dependencies/write_hif_targets/20200428.v1/HIF_Target_Genes.20200428.v1.tsv", data.table = F)
## input protein protein interaction table
ppi_df <- fread(data.table = F, input = "./Resources/Knowledge/PPI/protein_pair_table_v2.txt")
## input MYC targets
myc_tf_tab <- fread(input = "./Resources/Analysis_Results/dependencies/write_myc_targets/20200227.v1/MYC_Target_Genes.20200227.v1.tsv", data.table = F)
## input NRF2 targets
nrf2_tf_tab <- fread(input = "./Resources/Analysis_Results/dependencies/write_nrf2_targets/20200302.v1/NRF2_Target_Genes.20200302.v1.tsv", data.table = F)
## input TP53 targets
tp53_tf_tab <- fread(input = "./Resources/PPI/TF_interactions_TP53_manual.txt", data.table = F)
## input SQSTM1 (5q gain) affected genes
sqstm1_df <- data.frame(source_genesymbol = "SQSTM1",
                        target_genesymbol = c("KEAP1", "NFE2L2", 
                                              "TRAF6", "NFKB1", "MTOR"))

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
                                            filter(GENE %in% c("PBRM1", "SETD2", "BAP1", "KDM5C", "KDM6C")) %>%
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
genes_epigenetic_smgs_related_df$pathway_name <- "Epigenetic machinary"
# merge all pathways ------------------------------------------------------
genes_df <- rbind(genes_vhl_hif_df,
                  genes_epigenetic_smgs_related_df)

# write table -------------------------------------------------------------
write.table(x = genes_df, file = paste0(dir_out, "ccRCC_Pathogenic_Pathways_Genes.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
