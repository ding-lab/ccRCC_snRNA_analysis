# Yige Wu @ WashU 2021 Nov

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input ID mapping
idmapping_df <- readxl::read_xlsx(path = "~/Documents/Project/ccRCC_snRNA/Submission/CDS_upload/GSM_id_mapping.xlsx")

# process snATAC data -----------------------------------------------------
snATAC_metadata_df <- metadata_df %>%
  filter(Case != "C3L-00359") %>%
  filter(snATAC_used) %>%
  mutate(filename = paste0(Aliquot.snRNA, "_fragments.tsv.gz"))
idmapping_df = idmapping_df %>%
  mutate(aliquot_id = str_split_fixed(string = project_sample_id, pattern = "_", n = 2)[,1]) %>%
  filter(grepl(x = project_sample_id, "snATAC"))
  
snATAC_metadata_df$GEO_sample_id = mapvalues(x = snATAC_metadata_df$Aliquot.snRNA, from = idmapping_df$aliquot_id, to = as.vector(idmapping_df$GEO_sample_id))
final_df = snATAC_metadata_df %>%
  select(GEO_sample_id, filename)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "GEO.ccRCC.addfragments.", run_id, ".tsv")
write.table(x = final_df, file = file2write, quote = F, sep = "\t", row.names = F)


