# Yige Wu @WashU Sep 2019
## for summarizing xCell estimate of immune cell + stroma cell fractions from CPTAC ccRCC manuscript

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# Input xCell stats -------------------------------------------------------
xcell_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S7.xlsx", sheet = "xCell Signatures", skip = 2)
## because there is a row in character format so the rest of the data frame is in character format, need to change it double format
xcell_value_tab <- xcell_tab %>%
  filter(Samples != "Immune Group")

xcell_value_mat <- as.matrix(xcell_value_tab[,2:ncol(xcell_value_tab)]) %>% as.numeric()
rownames(xcell_value_mat) <- xcell_value_tab$Samples
xcell_value_mat %>%
  head()
# input table of cell family from xCell -----------------------------------



# summarize and write table -----------------------------------------------


