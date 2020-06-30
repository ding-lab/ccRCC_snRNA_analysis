# Yige Wu @WashU Feb 2020
## for writing a table for HIF targets

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
## input TF integrations
tf_tab <- fread(input = "./Resources/Knowledge/PPI/TF_interactions.txt", data.table = F)

# filter ------------------------------------
hif_tf_tab <- tf_tab %>%
  filter(source_genesymbol %in% c("HIF1A", "EPAS1")) %>%
  select(source_genesymbol, target_genesymbol) %>%
  mutate(target_genefunction = NA)

# add some canonical targets from literature ------------------------------
## ref1: https://www.nature.com/articles/nrc3183
hif_targets_manual1 <- readxl::read_excel(path = "./Resources/PPI/HIF_target_genes.xlsx", sheet = "Sheet1")
colnames(hif_targets_manual1)
hif_targets_manual1 <- hif_targets_manual1 %>%
  mutate(target_genesymbol = str_split_fixed(string = `Gene (protein)`, pattern = "\\s", n = 2)[,1])
hif_targets_manual12merge <- melt(hif_targets_manual1 %>%
                                    select(target_genesymbol, `HIF1α target gene`, `HIF2α target gene`), 
                                  id.vars = c("target_genesymbol"))
hif_targets_manual12merge <- hif_targets_manual12merge %>%
  filter(value == "+") %>%
  mutate(source_genesymbol = ifelse(variable == "HIF2α target gene", "EPAS1", "HIF1A")) %>%
  select(source_genesymbol, target_genesymbol)
## add Function
hif_targets_manual12merge$target_genefunction <- mapvalues(x = hif_targets_manual12merge$target_genesymbol, from = hif_targets_manual1$target_genesymbol, to = as.vector(hif_targets_manual1$Function))
## merge
hif_targets <- rbind(hif_tf_tab, hif_targets_manual12merge)

# add and correction gene functions ---------------------------------------
hif_targets$target_genefunction[hif_targets$target_genesymbol %in% c("VEGFA", "VEGFB", 
                                                                     "FLT1", "EGF", "SERPINE1", "ANGPT2", "TEK", "TIMP1",
                                                                     "PDGF", "ADM")] <- "Angiogenesis"
hif_targets$target_genefunction[hif_targets$target_genesymbol %in% c("TF", "TFRC")] <- "Iron metabolism"
hif_targets$target_genefunction[hif_targets$target_genesymbol %in% c("END1", "NOS2", "NOS3", "HMOX1", "NPPA")] <- "Vescular tone"
hif_targets$target_genefunction[hif_targets$target_genesymbol %in% c("LTBR", "TLR2")] <- "Inflammation"
hif_targets$target_genefunction[hif_targets$target_genesymbol %in% c( "HK1", "HK2", "HK3", "HKDC", "PFKL", 
                                                                      "GAPDH", "ALDOA", "ALDOB", "ALDOC",
                                                                      "ENO1", "ENO2", "ENO3", "ENO4", "PGK1", "PGK2",
                                                                      "PFKFB3", 
                                                                      "LDHA", "LDHB", "LDHC", "LDHAL6A", "LDHAL6B",
                                                                      "SLC2A1")] <- "Promote anaerobic metabolism"
hif_targets$target_genefunction[hif_targets$target_genesymbol %in% c( "PDK1")] <- "Inhibit TCA cycle metabolism"
hif_targets$target_genefunction[hif_targets$target_genesymbol %in% c("BNIP3", "BCL2", "RORA")] <- "Apoptosis Autophagy"
hif_targets$target_genefunction[hif_targets$target_genesymbol %in% c("CDKN1A", "CDKN1B")] <- "Cell cycle progression"
hif_targets$target_genefunction[hif_targets$target_genesymbol %in% c("MET", "EGFR", "TGFA", "TGFBR1", "ID2", "MYC")] <- "Proliferation survival"
## ID2: https://www.sciencedirect.com/science/article/pii/S2211383515000817
hif_targets$target_genefunction[hif_targets$target_genesymbol %in% c("CA9", "CA12")] <- "pH homeostasis"
hif_targets <- unique(hif_targets)

# write table -------------------------------------------------------------
write.table(x = hif_targets, file = paste0(dir_out, "HIF_Target_Genes.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

