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
## input known TF relationship
tf2target_df <- fread(data.table = F, input = "./Resources/Knowledge/PPI/Transcriptional/omnipathdb.transcriptional.20200908.txt")

# filter by DEGs and TFs ----------------------------------------------------------
tf2target_filtered_df <- tf2target_df %>%
  filter(source_genesymbol %in% c("RELA", "REL", "NFKB1", "NFKB2", "MXI1", "NFIC", "NFIA", "HSF2", "TP73", "NR3C2", "HIF1A")) %>%
  filter(target_genesymbol %in% c("PFKP", "GATM", "PVT1", "ERGIC1", "NDRG1", "LINC00887", "PLIN2", "PAM", 
                                  "TNIP1", "PHYKPL", "LRRC41", "EGLN3", "CDK18",
                                  "SEMA5B", "CIT", "ZNF395", "BIRC3", "C16orf74",
                                  "SPIRE1", "FRMD3", "NNMT", "FAM13A", "GIT2",
                                  "SEMA6A", "FOXP2", "COL23A1", "MXI1", "ST6GAL1",
                                  "EFNA5", "HAVCR1", "LUCAT1", "EVA1C", "KLF7",
                                  "PLOD2", "PLSCR1", "CDH6", "DNAH11", "PDK1",
                                  "PKM", "ZNF608", "GRAMD2B", "TNFAIP8", "LINCO1426",
                                  "SLC38A1", "CP", "SLC6A3", "SGPP2", "ENPP3", 
                                  "MYOCOS", "KCTD3", "TSC22D3", "BARX2", "GAS2L3",
                                  "MSC-AS1", "PTGER3", "SLC16A1-AS1", "HILPDA", "STK39",
                                  "EGFR", "PCSK6", "RETREG1", "LINC02471", "TRIM9",
                                  "TMEM232", "SPAG17", "ABI3BP", "NREP", "RNF145", 
                                  "PRELID2", "ABLIM3", "APBB11P", "MGAM", "TSPAN12", 
                                  "PLCB1")) %>%
  arrange(target_genesymbol)


