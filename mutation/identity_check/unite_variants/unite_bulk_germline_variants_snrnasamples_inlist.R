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
dir_vcf <- "./Resources/Bulk_Processed_Data/Germline_Variants/VCF_Files_pickCaller20200706/"
filenames_vcf <- list.files(path = dir_vcf)
filenames_vcf <- filenames_vcf[grepl(x = filenames_vcf, pattern = "vcf")]
filenames_vcf
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## get the cases to process
ids_case2process <- unique(idmetadata_df$Case[idmetadata_df$snRNA_available])

# input vcf file and process ----------------------------------------------
vcf_bulk_list <- list()
for (filename_vcf in filenames_vcf) {
  ## get aliquot id
  id_case <- str_split_fixed(string = filename_vcf, pattern = "\\.", n = 3)[,1]
  if (!(id_case %in% ids_case2process)) {
    next()
  }
  ## get the path to the VCF file
  path_vcf <- paste0(dir_vcf, filename_vcf)
  vcf_df <- fread(data.table = F, skip = "#CHROM", input = path_vcf)
  if (nrow(vcf_df) == 0) {
    next()
  }

  ## rename last column and add aliquot id
  vcf_df <- vcf_df %>%
    rename('CHROM' = '#CHROM') %>%
    mutate(ID_VARIANT = paste0(CHROM, "_", POS, "_", REF, "_", ALT))
  vcf_bulk_list[[id_case]] <- vcf_df$ID_VARIANT
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Germline_Variants.", "bulk.", "pickCaller.", "List.", run_id, ".RDS")
saveRDS(object = vcf_bulk_list, file = file2write, compress = T)
