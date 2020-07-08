# Yige Wu @WashU Jul 2020

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

# input dependencies -----------------------------------------------------
## input VAFs for bulk data
vaf_bulk_df <- fread(data.table = F, input = "./Resources/Analysis_Results/mutation/identity_check/calculate_vaf/calculate_vaf_for_bulk_germline_variants/20200706.v1/Germline_Variants.bulk.pickCaller.Shared_Variants.snATAC_Cases.VAF.20200706.v1.tsv")
## input VAFs for snATAC data
vaf_snatac_df <- fread(data.table = F, input = "./Resources/Analysis_Results/mutation/identity_check/calculate_vaf/calculate_vaf_for_snatac_germline_variants/20200706.v1/Germline_Variants.snATAC.pickCaller.Shared_Variants.VAF.20200706.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")

# identify shared variants between bulk and snatac------------------------------------------------
## filter by VAF
vaf_bulk_df <- vaf_bulk_df %>%
  filter(!is.na(VAF))
vaf_snatac_df <- vaf_snatac_df %>%
  filter(!is.na(VAF))
variants_bulk <- unique(vaf_bulk_df$ID_VARIANT)
length(variants_bulk)
variants_snatac <- unique(vaf_snatac_df$ID_VARIANT)
length(variants_snatac)
variants_shared <- intersect(variants_bulk, variants_snatac)
length(variants_shared)
variants_shared %>% head()

# make data frame for calculating correlation -----------------------------
## add readable aliquot id to the vaf_snatac_df
vaf_snatac_filtered_df <- vaf_snatac_df %>%
  filter(ID_VARIANT %in% variants_shared) %>%
  arrange(ID_VARIANT)
vaf_snatac_filtered_df$Id_Aliquot_WU <- mapvalues(x = vaf_snatac_filtered_df$ID_ALIQUOT, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
vaf_snatac_filtered_df <- vaf_snatac_filtered_df %>%
  mutate(Id_Aliquot_Data = paste0("snATAC.", Id_Aliquot_WU))
unique(vaf_snatac_filtered_df$Id_Aliquot_WU)
vaf_bulk_filtered_df <- vaf_bulk_df %>%
  filter(ID_VARIANT %in% variants_shared) %>%
  mutate(Id_Aliquot_Data = paste0("bulk.", ID_Case, "-T1")) %>%
  arrange(ID_VARIANT)
## combine
vaf_combined_df <- rbind(vaf_snatac_filtered_df %>%
                           select(ID_VARIANT, Id_Aliquot_Data, VAF), 
                         vaf_bulk_filtered_df%>%
                           select(ID_VARIANT, Id_Aliquot_Data, VAF))

vaf_combined_dcast_df <- dcast(data = vaf_combined_df,
                                    formula = ID_VARIANT ~ Id_Aliquot_Data, value.var = c("VAF"))
nrow(vaf_combined_dcast_df)
## remove rows with any NAs
vaf_combined_dcast_df <- vaf_combined_dcast_df[rowMeans(is.na(vaf_combined_dcast_df[,-1])) == 0,]
nrow(vaf_combined_dcast_df)
## 2334 rows
# get correlation result --------------------------------------------------
cor_coef <- cor(vaf_combined_dcast_df[,-1])
cor_coef %>% head()
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Correlation_Coeffcients.", "snATAC_vs_bulk.", "snATAC_Cases.", run_id, ".RDS")
saveRDS(object = cor_coef, file = file2write, compress = F)

