# set shared snRNA processing parameters ----------------------------------
num_pc <- 30
num_features <- 3000

# gene lists --------------------------------------------------------------
## significantly mutated genes
ccRCC_SMGs <- c("VHL", 
                "PBRM1", 
                "SETD2", "BAP1", 
                "KDM5C", "KDM6A",
                "PTEN", 
                "MTOR", "TSC1",
                "TP53")
## PBAF gens
### reference: https://www.nature.com/articles/onc20094/figures/1
### reference: https://www.nature.com/articles/onc20094/tables/1
pbaf_genes <- c("PBRM1", "ARID2", "SMARCE1", "SMARCC2", "ACTL6A", "ACTL6B", "SMARCC1", "SMARCD1", "SMARCB1")

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
                "SQSTM1",
                "FGFR4")
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
                "5q", "5q")
ccrcc_cna_genes_df <- data.frame(gene_symbol = c(genes_loss, genes_gain),
                                 gene_cna_type = c(rep("Loss", length(genes_loss)), rep("Gain", length(genes_gain))),
                                 chr_region = chr_region)
ccrcc_cna_genes_df
