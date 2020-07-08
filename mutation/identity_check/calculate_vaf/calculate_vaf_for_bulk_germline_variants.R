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
## input bulk germline variants
vcf_df <- fread(data.table = F, input = "./Resources/Analysis_Results/mutation/identity_check/intersect_bulk_germline_variants_for_snatac_samples/20200706.v1/Germline_Variants.bulk.pickCaller.Shared_Variants.snATAC_Cases.20200706.v1.tsv")
length(unique(vcf_df$ID_Case))
length(unique(vcf_df$ID_VARIANT))

# get allelic depths for the alt alleles and total depth --------------------------
## get alterative allele count
vcf_df$Alt_Count <- sapply(X = 1:nrow(vcf_df), FUN = function(i, vcf_data_frame) {
  string_format <- vcf_data_frame[i, "FORMAT"]
  string_format_split <- str_split(string = string_format, pattern = ":")[[1]]
  string_genotype <- vcf_data_frame[i, "GENOTYPE"]
  string_genotype_split <- str_split(string = string_genotype, pattern = ":")[[1]]
  if (length(which(string_format_split == "AD")) == 0) {
    result <- NA
  } else {
    string_AD <- string_genotype_split[which(string_format_split == "AD")]
    string_AD_split <- str_split(string = string_AD, pattern = ",")[[1]]
    result <- string_AD_split[length(string_AD_split)]
  }
  return(result)
}, vcf_data_frame = vcf_df)
vcf_df$Alt_Count <- as.numeric(vcf_df$Alt_Count)
vcf_df$Alt_Count
## get total allele count
vcf_df$Total_Count <- sapply(X = 1:nrow(vcf_df), FUN = function(i, vcf_data_frame) {
  string_format <- vcf_data_frame[i, "FORMAT"]
  string_format_split <- str_split(string = string_format, pattern = ":")[[1]]
  string_genotype <- vcf_data_frame[i, "GENOTYPE"]
  string_genotype_split <- str_split(string = string_genotype, pattern = ":")[[1]]
  if (length(which(string_format_split == "DP")) == 0) {
    result <- NA
  } else {
    string_DP <- string_genotype_split[which(string_format_split == "DP")]
    result <- string_DP
  }
  return(result)
}, vcf_data_frame = vcf_df)
vcf_df$Total_Count <- as.numeric(vcf_df$Total_Count)
vcf_df$Total_Count
## calculate VAF
vcf_df <- vcf_df %>%
  mutate(VAF = Alt_Count/Total_Count)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Germline_Variants.", "bulk.", "pickCaller.", "Shared_Variants.", "snATAC_Cases.", "VAF.", run_id, ".tsv")
write.table(x = vcf_df, file = file2write, quote = F, sep = "\t", row.names = F)




