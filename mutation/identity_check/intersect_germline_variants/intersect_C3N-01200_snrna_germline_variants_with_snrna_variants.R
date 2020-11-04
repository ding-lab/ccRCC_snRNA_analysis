# Yige Wu @WashU OCt 2020

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
vcf_united_list <- readRDS(file = "./Resources/Analysis_Results/mutation/identity_check/unite_variants/unite_snrna_germline_variants_all_samples_inlist/20201029.v1/Germline_Variants.snRNA.pickCaller.List.20201029.v1.RDS")

# make loop ---------------------------------------------------------------
aliquot1_vec <- NULL
aliquot2_vec <- NULL
number_variants_vec <- NULL
i <- 0
for (aliquot1 in c("C3N-01200-T1", "C3N-01200-T2", "C3N-01200-T3", "C3N-01200-N")) {
  for (aliquot2 in names(vcf_united_list)) {
    variants_intersect_tmp <- intersect(vcf_united_list[[aliquot1]], vcf_united_list[[aliquot2]])
    number_variants_tmp <- length(variants_intersect_tmp)
    ## store results
    i <- i+1
    aliquot1_vec[i] <- aliquot1
    aliquot2_vec[i] <- aliquot2
    number_variants_vec[i] <- number_variants_tmp
  }
}
## make data frame
number_variants_df <- data.frame(aliquot1 = aliquot1_vec,
                                 aliquot2 = aliquot2_vec,
                                 number_variants = number_variants_vec)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Number_Variants.C3N-01200_snRNA_vs_All_snRNA", ".tsv")
write.table(x = number_variants_df, file = file2write, quote = F, sep = "\t", row.names = F)
