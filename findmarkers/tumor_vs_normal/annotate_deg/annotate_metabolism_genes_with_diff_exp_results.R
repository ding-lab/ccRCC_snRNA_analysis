# Yige Wu @WashU Jun 2021
## gene reference: TCGA ccRCC study Fig. 5b

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
## input united DEG result
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_bulkRNA_protein_DEGs/20210601.v1/Tumor_vs_PT_DEGs.United.snRNA.bulkRNA.Protein.20210601.v1.tsv")
## input the metabolism genes to look at
genes_process <- c("HK1", "HK2", "HK3", "GPI", "PFKP", "PFKL", "PFKM", "ALDOB", "ALDOA", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGK2",
                   "PGAM1", "PGAM2", "ENO1", "ENO2", "ENO3", "PKM", "PKLR", "LDHA", "LDHB", "LDHC", "LDHD",
                   "PDK1", "PDK2", "PDK3", "PDK4", "PDP1", "PDP2", "DLAT", "DLD", "PDHA1", "PDHA2", "PDHB",
                   "ACO2", "IDH2", "OGDH", "SDHB", "SDHC", "SDHD", "FH", "MDH2",
                   "ACLY", "ACACA", "FASN",
                   "PGM1", "PGM2", "UGP2", "GYS1", "PYGL", "GBE1")


# filter ------------------------------------------------------------------
degs_filtered_df <- degs_df %>%
  filter(genesymbol_deg %in% genes_process)
rownames(degs_filtered_df) <- degs_filtered_df$genesymbol_deg
rownames_ordered <- genes_process[genes_process %in% degs_filtered_df$genesymbol_deg]
degs_filtered_df <- degs_filtered_df[rownames_ordered,]
degs_filtered_selected_df <- degs_filtered_df %>%
  select(genesymbol_deg, Num_sig_up, Num_sig_down, logFC.bulkRNA, FDR.bulkRNA, direction.bulkRNA, meddiff_exp.bulkpro, FDR.bulkpro, direction.bulkpro)
