# Yige Wu @ WashU 2020 Dec

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/load_data.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input mutation table from the data freeze
maf_raw_df <- fread(data.table = F, input = "../ccRCC_Drug/Resources/Bulk_Processed_Data/Data_Files/batch8_cellline/somaticMut/somaticMut_cellLine_hk2.dnp.meta2.non-silent.tsv")
## function
get_somatic_mutation_aachange_vaf_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat & maf$Variant_Classification != "Silent",]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  maf$vaf <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
  maf$HGVSp_sim <- gsub(x = maf$HGVSp_Short, pattern = "p\\.", replacement = "")
  
  maf$aachange_vaf <- paste0(maf$HGVSp_sim, "(", signif(x = maf$vaf, digits = 2), ")")
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    aahange_vaf <- paste0(unique(x), collapse = ",")
    return(aahange_vaf)
  }, value.var = "aachange_vaf", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

# make matrix for aa change & vaf------------------------------------------------------------
mut_mat <- get_somatic_mutation_aachange_vaf_matrix(pair_tab = ccRCC_SMGs, maf = maf_raw_df)
mut_df <- t(mut_mat[,-1]) %>% as.data.frame()
## order
genesymbols_mut <- colnames(mut_df)
mut_df$Case <- rownames(mut_df)
genesymbols_mut_driver <- ccRCC_drivers[ccRCC_drivers %in% genesymbols_mut]
genesymbols_mut_othersmg <- genesymbols_mut[!(genesymbols_mut %in% genesymbols_mut_driver)]
mut_df <- mut_df[, c("Case", genesymbols_mut_driver, genesymbols_mut_othersmg)]

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "snRNA_ccRCC_Mutation_AAchange_VAF_Table.", run_id, ".tsv")
write.table(x = mut_df, file = file2write, sep = "\t", row.names = F, col.names = T, quote = F)
