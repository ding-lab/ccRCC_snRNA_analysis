# Yige Wu @WashU March 2020
## running on local
## for merging the bulk omics profile and snRNA-based CNV
## bulk omics profile include all the  mutation (SMG), CNV (reported by TCGA, bulk), VHL methylation and 3p translocation

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
## input bulk omics profile
## input te bulk genomics/methylation events
merged_bulk_events_df <- fread(input = "./Resources/Analysis_Results/bulk/other/merge_bulk_events/20210504.v1/merged_bulk_events.20210504.v1.tsv", data.table = F)
## input snRNA copy number profile for tumor cells
merged_sn_events_df <- fread(input = "./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/estimate_fraction_of_tumorcells_with_expectedcnv_perchrregion_per_sample_using_cnvgenes/20200318.v1/fraction_of_tumorcells.expectedCNA.by_chr_region.20200318.v1.tsv", data.table = F)
## input the tumor content
merged_tumorpurity_df <- fread(input = "./Resources/Analysis_Results/bulk/tumor_content/merge_tumor_content_from_bulk_and_snRNA/20200319.v1/Perc_Tumor_Content_from_snRNA_and_bulkRNA.20200319.v1.tsv", data.table = F)
## input id meta data
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv", data.table = F)

# filter snRNA copy number to 3p, 5q, 14q ---------------------------------
merged_sn_events_df <- merged_sn_events_df %>%
  rename(CN.sn.3p_loss.fraction = `3p`) %>%
  rename(CN.sn.5q_gain.fraction = `5q`) %>%
  rename(CN.sn.14q_loss.fraction = `14q`) %>%
  select(aliquot, CN.sn.3p_loss.fraction, CN.sn.5q_gain.fraction, CN.sn.14q_loss.fraction)
## add case ids for the aliquot ids
merged_sn_events_df$Case <- mapvalues(x = merged_sn_events_df$aliquot, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))
merged_sn_events_df$Aliquot_snRNA_WU <- mapvalues(x = merged_sn_events_df$aliquot, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Aliquot.snRNA.WU))

# merge sn cnv with bulk omics profile ------------------------------------
## filter meta data by just the original tumor segment (discovery cohort)
id_metadata_original_df <- id_metadata_df %>%
  filter(Is_discovery_set == T) %>%
  filter(Sample_Type == "Tumor")
## merge by snRNA aliquot ids
bulk_sn_omicsprofile_df <- merge(merged_bulk_events_df, merged_sn_events_df, 
                                 by.x = c("Case", "Aliquot_snRNA_WU"), by.y = c("Case", "Aliquot_snRNA_WU"), all = T)
bulk_sn_omicsprofile_df <- bulk_sn_omicsprofile_df %>%
  rename(Aliquot.snRNA = aliquot)
## merge with bulk and sn tumor purity
bulk_sn_omicsprofile_df <- merge(bulk_sn_omicsprofile_df, 
                                 merged_tumorpurity_df %>%
                                   select(Case, Aliquot, Frac_Tumor_Barcodes, ESTIMATE_TumorPurity_RNA) %>%
                                   rename(Aliquot.snRNA = Aliquot) %>%
                                   rename(TumorPurity.snRNA = Frac_Tumor_Barcodes) %>%
                                   rename(TumorPurity.bulk = ESTIMATE_TumorPurity_RNA), 
                                 by = c("Case", "Aliquot.snRNA"), 
                                 all = T)
bulk_sn_omicsprofile_df$Aliquot.snRNA <- mapvalues(x = bulk_sn_omicsprofile_df$Aliquot_snRNA_WU, from = id_metadata_df$Aliquot.snRNA.WU, to = as.vector(id_metadata_df$Aliquot.snRNA))

# write table -------------------------------------------------------------
file2write <- paste0(dir_out, "bulk_sn_omics_profile.", run_id, ".tsv")
write.table(x = bulk_sn_omicsprofile_df, file = file2write, quote = F, sep = "\t", row.names = F)

