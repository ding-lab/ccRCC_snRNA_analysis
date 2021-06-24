# Yige Wu @WashU Dec 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/load_data.R")

## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the bulk mutation data
maf_df <- loadMaf()
# maf_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Somatic_mutation/WashU/CPTAC_ccRCC_discovery_somatic_mutation.maf_v1.0.tsv")
## input BAP1 related proteins
bap1_interactors_df <- fread(data.table = F, input = "../ccRCC_confirmatory/Resources/Database/Protein_Protein_Interactions/BioGRID/BIOGRID-GENE-113911-4.3.196.tab3.txt")

# annotate samples based on PBRM1 & BAP1 mutation -------------------------
maf_case_df <- maf_df %>%
  filter(Tumor_Sample_Barcode == "C3L-00096_T")

maf_case_bap1interactor_df <- maf_case_df %>%
  filter(Hugo_Symbol %in% c(bap1_interactors_df$`Official Symbol Interactor A`, bap1_interactors_df$`Official Symbol Interactor B`))


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "C3L-00096_T.Mutations.", run_id, ".tsv")
write.table(x = maf_case_df, file = file2write, quote = F, sep = "\t", row.names = F)

