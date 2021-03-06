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
dir_vcf <- "./Resources/snRNA_Processed_Data/Germline_Variants_snRNA/VCF_Files_pickCaller20200706/"
filenames_vcf <- list.files(path = dir_vcf)
filenames_vcf <- filenames_vcf[grepl(x = filenames_vcf, pattern = "vcf")]
filenames_vcf
## input snRNA samples
snatacsample_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Data_Generation/snATAC_Sample_Selection_Rationale.tsv")
# input vcf file and process ----------------------------------------------
vcf_merged_df <- NULL
variants_intersect <- NULL
for (filename_vcf in filenames_vcf) {
  ## get aliquot id
  ## get aliquot id
  id_aliquot <- str_split_fixed(string = filename_vcf, pattern = "\\.", n = 3)[,1]
  id_aliquot
  if (!(id_aliquot %in% snatacsample_df$Aliquot)) {
    next()
  }
  path_vcf <- paste0(dir_vcf, filename_vcf)
  vcf_df <- fread(data.table = F, skip = "#CHROM", input = path_vcf)
  if (nrow(vcf_df) == 0) {
    next()
  }

  ## rename last column and add aliquot id
  vcf_df <- vcf_df %>%
    rename("GENOTYPE" = id_aliquot) %>%
    rename('CHROM' = '#CHROM') %>%
    mutate(ID_ALIQUOT = id_aliquot) %>%
    mutate(ID_VARIANT = paste0("chr", CHROM, "_", POS, "_", REF, "_", ALT))
  vcf_merged_df <- rbind(vcf_merged_df, vcf_df)
  ## merge variants
  if (!is.null(variants_intersect)) {
    variants_intersect <- intersect(variants_intersect, vcf_df$ID_VARIANT)
  } else {
    variants_intersect <- vcf_df$ID_VARIANT
  }
}
length(variants_intersect)

# filter down to only shared variants -------------------------------------
vcf_shared_df <- vcf_merged_df %>%
  filter(ID_VARIANT %in% variants_intersect)
nrow(vcf_shared_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Germline_Variants.", "snRNA.", "pickCaller.", "Shared_Variants.", "snATAC_Samples.", run_id, ".tsv")
write.table(x = vcf_shared_df, file = file2write, quote = F, sep = "\t", row.names = F)
