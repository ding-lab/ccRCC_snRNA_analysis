# set shared snRNA processing parameters ----------------------------------
num_pc <- 30
num_features <- 3000

# gene lists --------------------------------------------------------------
## SMG source 1: Sato et al, Nature Genetics, 2013
## SMG source 1: TCGA, Nature, 2013
ccRCC_SMGs <- c("VHL", "PBRM1", "BAP1", "TCEB1", "SETD2", "TP53", "PIK3CA", "MTOR", "KDM5C", "PTEN","TSC1",
                "CCNB2",
                "MSR1", "TXNIP", "BTNL3", "SLITRK6", "RHEB", "ARID1A", "NPNT", 
                "FPGT", "MUDENG", "TET2", "MUC4", "MLLT10", "MSGN1", "KRT32", "M6PR", "RPL14", "GRB7", "CSMD3", "DNHD1", "NLRP12", "VMO1", "OR4C13", "KCNMA1", "LMAN2L", "ZNF536", "YIPF3")
# Two examples
# of cullin–RING ubiquitin ligase system molecular assemblies using CUL2 or CUL5 (left) and CUL3 (right) that interact with the BC-box protein–
# Elongin C–Elongin B complex and BTB protein, respectively, to recruit substrate for ubiquitination and subsequent degradation. VHL and KEAP1 are examples of BC-box and BTB proteins, respectively, that recruit HIF and NRF2 proteins for ubiquitin-mediated degradation. 
vhl_complex_genes <- c("VHL", 
                       "TCEB1", "TCEB2",
                       "ELOC", "ELOB",
                       "RBX1", "RNF7",
                       "CUL2", "CUL5")
keap1_nrf2_complex_genes <- c("NFE2L2", "KEAP1", "CUL3", "RBX1")
ubiquitin_proteasome_genes <- c(vhl_complex_genes,
                                keap1_nrf2_complex_genes,
                                "USP24", "NEDD4", "WWP2", "URB1", "USP34")
swisnf_genes <- c("PBRM1", "ARID1A", "ARID1B", "ARID2", "SMARCA2", "SMARCA4", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD2", "SMARCD2")
other_epigeneticregulator_genes <- c("SETD2", 
                                     "BAP1", 
                                     "KDM5C", "KDM6A",
                                     "TET2")
pi3k_mtor_genes <- c("EGFR", "ERBB3", "FGFR3", "FGFR4", "IGF1R", "PIK3CA", "PIK3CB", "PIK3CG","PTEN", "AKT1", "AKT2", "AKT3", "MTOR", "TSC1", "TSC2", "RPS6KA2", "RPS6KA3", "RPS6KA6", "RHEB", 
                     "SRC", "PKT2")
p53_cellcycle_genes <- c("ATM", "CHECK2", "TP53", "MDM2", "CDKN2A", "CCND1", "E2F3", "CCNB2",
                         "MYC")
fat_cadherins_genes <- c("FAT1", "FAT2", "FAT3", "FAT4")

## genes mutated in ccRCC
genes_pathogenicpathways_in_ccRCC <- c(ubiquitin_proteasome_genes,
                            swisnf_genes,
                            other_epigeneticregulator_genes,
                            pi3k_mtor_genes,
                            p53_cellcycle_genes,
                            fat_cadherins_genes)
genes_pathogenicpathways_in_ccRCC <- unique(genes_pathogenicpathways_in_ccRCC)
## PBAF gens
### reference: https://www.nature.com/articles/onc20094/figures/1
### reference: https://www.nature.com/articles/onc20094/tables/1
pbaf_genes <- c("PBRM1", "ARID2", "SMARCE1", "SMARCC2", "ACTL6A", "ACTL6B", "SMARCC1", "SMARCD1", "SMARCB1")


# make known CNV regions --------------------------------------------------
genes_loss <- c("VHL", "PBRM1", "BAP1", "SETD2",
                "HIF1A",
                "CDKN2A", 
                "PTEN", 
                "NEGR1",
                "QKI",
                "CADM2", 
                "PTPRD", 
                "NRXN3")
genes_gain <- c("PRKCI", "MECOM",
                "MDM4",
                "MYC",
                "JAK2",
                "SQSTM1","FGFR4", "EZH2", "STC2", "VCAN") # https://www.sciencedirect.com/science/article/pii/S0304419X18300441?casa_token=E86Vf15p0XAAAAAA:xFno1e29kPmGAtLoKWPmqHktdYNj5OEGO9Klz1svFjrgH5DkqsLWGjlDkynI4LM4998BFN4D
chr_region <- c(rep("3p", 4),
                "14q",
                "9p21", 
                "10q23",
                "1p31",
                "6q24",
                "3p12",
                "9p23",
                "14q24",
                "3p26", "3p26",
                "1q32",
                "8q24",
                "9q24",
                rep("5q", 5))
ccrcc_cna_genes_df <- data.frame(gene_symbol = c(genes_loss, genes_gain),
                                 gene_cna_type = c(rep("Loss", length(genes_loss)), rep("Gain", length(genes_gain))),
                                 chr_region = chr_region)
ccrcc_cna_genes_df
