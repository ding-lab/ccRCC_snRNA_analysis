# Yige Wu @WashU Nov 2020
## 2021-08-05 added two more tumor samples

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
## input id meta daa
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")

# make paths --------------------------------------------------------------
dir_base_katmai <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
srat_paths_df <- idmetadata_df %>%
  filter(snRNA_available) %>%
  filter(Sample_Type == "Tumor") %>%
  select(Case, Sample, Sample_Type, Aliquot.snRNA, Aliquot.snRNA.WU) %>%
  mutate(Version = ifelse(Aliquot.snRNA.WU %in% c("C3N-00437-T1", "C3N-00317-T1"), "20210505.v1", "20201124.v1")) %>%
  mutate(Path_katmai = paste0(dir_base_katmai, "Data_Freezes/V2/snRNA/Tumor_Cell_Reclustered/", Aliquot.snRNA.WU, ".tumorcellreclustered.", 
                              Version, ".RDS"))
  
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Paths_TumorCellOnlyReclustered_SeuratObject.", run_id, ".tsv")
write.table(x = srat_paths_df, file = file2write, quote = F, sep = "\t", row.names = F)
