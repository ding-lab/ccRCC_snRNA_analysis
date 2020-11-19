# Yige Wu @WashU Nov 2020

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
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input samples whose cell type has been corrected
aliquots_reprocessed_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Individual_AllClusters/Cells_BySampleByClusterByCellTypeShorter.Over50.20201027.xlsx", sheet = "Sheet1")

# get aliquot ids --------------------------------------------------------------
## path for the 20 corrected tumor samples reclustered 20201109.v1
aliquots_reprocessed <- unique(aliquots_reprocessed_df$aliquot[aliquots_reprocessed_df$Cell_type.shorter.original == "Tumor cells" | aliquots_reprocessed_df$Cell_group.detailed == "Tumor cells"])
aliquots_reprocessed
## path for the  9 uncanged tumor cells reclustered version. 20200225.v1
idmetadata_filtered_df <- idmetadata_df %>%
  filter(Sample_Type == "Tumor") %>%
  filter(snRNA_available)
aliquots_unchanged <- idmetadata_filtered_df$Aliquot.snRNA[!(idmetadata_filtered_df$Aliquot.snRNA %in% aliquots_reprocessed)]
aliquots_unchanged

# make paths --------------------------------------------------------------
dir_base_katmai <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
srat_paths_df <- idmetadata_filtered_df %>%
  select(Case, Sample, Sample_Type, Aliquot.snRNA, Aliquot.snRNA.WU) %>%
  mutate(Path_katmai = ifelse(Aliquot.snRNA %in% aliquots_reprocessed, 
                              paste0(dir_base_katmai, "Data_Freezes/V1/snRNA/Tumor_Cell_Reclustered/", Aliquot.snRNA, ".malignant_nephron_epithelium.20201109.v1.RDS"),
                              paste0(dir_base_katmai, "Data_Freezes/V1/snRNA/Tumor_Cell_Reclustered/", Aliquot.snRNA, ".malignant_nephron_epithelium.20200225.v1.RDS"))) %>%
  mutate(Version = ifelse(Aliquot.snRNA %in% aliquots_reprocessed, "20201109.v1", "20200225.v1"))
                              
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Paths_TumorCellOnlyReclustered_SeuratObject.", run_id, ".tsv")
write.table(x = srat_paths_df, file = file2write, quote = F, sep = "\t", row.names = F)
