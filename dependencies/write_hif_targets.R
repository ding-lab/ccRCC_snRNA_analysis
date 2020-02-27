# Yige Wu @WashU Feb 2020
## for writing a table for HIF targets

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input TF integrations
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/PPI/TF_interactions.txt", data.table = F)

# filter ------------------------------------
hif_tf_tab <- tf_tab %>%
  filter(source_genesymbol %in% c("HIF1A", "EPAS1")) %>%
  select(source_genesymbol, target_genesymbol)

# add some canonical targets from literature ------------------------------
## ref1: https://www.nature.com/articles/nrc3183
hif_targets_manual1 <- readxl::read_excel(path = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/PPI/HIF_target_genes.xlsx", sheet = "Sheet1")
colnames(hif_targets_manual1)
hif_targets_manual1 <- hif_targets_manual1 %>%
  mutate(target_genesymbol = str_split_fixed(string = `Gene (protein)`, pattern = "\\s", n = 2)[,1]) %>%
  select(target_genesymbol, `HIF1α target gene`, `HIF2α target gene`)
hif_targets_manual12merge <- melt(hif_targets_manual1, id.vars = c("target_genesymbol"))
hif_targets_manual12merge <- hif_targets_manual12merge %>%
  filter(value == "+") %>%
  mutate(source_genesymbol = ifelse(variable == "HIF2α target gene", "EPAS1", "HIF1A")) %>%
  select(source_genesymbol, target_genesymbol)
## merge
hif_targets <- rbind(hif_tf_tab, hif_targets_manual12merge)
hif_targets <- unique(hif_targets)
# write table -------------------------------------------------------------
write.table(x = hif_targets, file = paste0(dir_out, "HIF_Target_Genes.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

