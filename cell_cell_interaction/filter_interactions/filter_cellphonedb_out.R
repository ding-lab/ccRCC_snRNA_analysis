# Yige Wu @WashU Sep 2020
## reference of the cellphonedb output: https://www.cellphonedb.org/documentation

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(monocle)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cellphonedb output
cellphone_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/CellPhoneDB/cell.phone.res.total.run20200818.txt")

# filter by cell count and filter out unknown----------------------------------------------------
cellphone_filtered_df <- cellphone_df %>%
  filter(cell.num.ct1 >= 10 & cell.num.ct2 >= 10) %>%
  filter(Cell_type1 != "Unknown" & Cell_type2 != "Unknown")

# filter out duplicates ---------------------------------------------------
cellphone_filtered_df <- cellphone_filtered_df %>%
  filter(!(interacting_pair %in% c("FLT1 complex_VEGFA", "FLT1 complex_VEGFB")))
  

# filter out cell-cell interactions that does not make sense --------------
# FGFR2_EPHA4
## A yeast two-hybrid analysis has shown that the juxtamembrane region of FGF receptor 3 (FGFR3) interacts with the cytoplasmic domain of EphA4, which is a member of the largest family of receptor tyrosine kinases. Complex formation between the two receptors was shown to be mediated by direct interactions between the juxtamembrane domain of FGFR1, FGFR2, FGFR3, or FGFR4 and the N-terminal portion of the tyrosine kinase domain of EphA4. 
### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1323220/
# LGALS9_MET
## cannot find the interaction in the innateDB
# VEGFA_GRIN2B
## cannot find any literature supporting this nor the I2D database
# AXL_IL15RA
## the original paper is retrated: https://pubmed.ncbi.nlm.nih.gov/16308569/
cellphone_filtered_df <- cellphone_filtered_df %>%
  filter(!(interacting_pair %in% c("FGFR2_EPHA4", "FGFR3_EPHA4", "FGFR4_EPHA4", "FGFR1_FGFR2") & as.vector(Cell_type1) != as.vector(Cell_type2))) %>%
  filter(!(interacting_pair %in% c("LGALS9_MET", "VEGFA_GRIN2B", "AXL_IL15RA")))


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cell.phone.res.total.run20200818.filtered.txt")
write.table(x = cellphone_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)