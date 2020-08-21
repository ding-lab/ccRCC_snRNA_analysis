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
## input VAFs for snRNA data
vaf_snrna_df <- fread(data.table = F, input = "./Resources/Analysis_Results/mutation/identity_check/calculate_vaf/calculate_vaf_for_snrna_germline_variants_shared_among_selectedsamples/20200820.v2/Germline_Variants.snRNA.pickCaller.Shared_Variants.SelectedSamples.VAF.20200820.v2.tsv")

# make data frame for calculating correlation -----------------------------
## add readable aliquot id to the vaf_snrna_df
vaf_snrna_filtered_df <- vaf_snrna_df
vaf_snrna_filtered_df$Id_Aliquot_WU <- mapvalues(x = vaf_snrna_filtered_df$ID_ALIQUOT, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
vaf_snrna_filtered_df <- vaf_snrna_filtered_df %>%
  mutate(Id_Aliquot_Data = paste0("snRNA.", Id_Aliquot_WU))
unique(vaf_snrna_filtered_df$Id_Aliquot_WU)
## combine
vaf_combined_df <- vaf_snrna_filtered_df
vaf_combined_dcast_df <- dcast(data = vaf_combined_df,
                                    formula = ID_VARIANT ~ Id_Aliquot_Data, value.var = c("VAF"))
nrow(vaf_combined_dcast_df)
## remove rows with any NAs
vaf_combined_dcast_df <- vaf_combined_dcast_df[rowMeans(is.na(vaf_combined_dcast_df[,-1])) == 0,]
nrow(vaf_combined_dcast_df)

# get correlation result --------------------------------------------------
cor_coef <- cor(vaf_combined_dcast_df[,-1])
cor_coef %>% head()
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Correlation_Coeffcients.", "snRNA.", "SelectedSamples.", run_id, ".RDS")
saveRDS(object = cor_coef, file = file2write, compress = F)

