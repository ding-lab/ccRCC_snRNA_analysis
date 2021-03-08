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
## input maf file
maf_df <- loadMaf()
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## function
get_somatic_mutation_aachange_vaf_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat & maf$Variant_Classification != "Silent",]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  maf$vaf <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
  
  maf$aachange_vaf <- paste0(maf$HGVSp_Short, "(", signif(x = maf$vaf, digits = 2), ")")
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    aahange_vaf <- paste0(unique(x), collapse = ",")
    return(aahange_vaf)
  }, value.var = "aachange_vaf", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

# set parameters ----------------------------------------------------------
## set genes to look for mutations
genes2process <- ccRCC_SMGs
## set sample case ids to process
caseids2process <- unique(idmetadata_df$Case[idmetadata_df$Aliquot.snRNA.WU %in% snatacsample_anno_df$easyid])
## filter maf file
maf_filtered_df <- maf_df %>%
  mutate(id_case = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", 2)[,1]) %>%
  filter(id_case %in% caseids2process)

# make matrix for vaf-------------------------------------------------------------
vaf_mat <- get_somatic_mutation_vaf_matrix(pair_tab = genes2process, maf = maf_filtered_df)
vaf_df <- t(vaf_mat[,-1]) %>% as.data.frame()
## order
genesymbols_mut <- colnames(vaf_df)
vaf_df$Case <- rownames(vaf_df)
genesymbols_mut_driver <- ccRCC_drivers[ccRCC_drivers %in% genesymbols_mut]
genesymbols_mut_othersmg <- genesymbols_mut[!(genesymbols_mut %in% genesymbols_mut_driver)]
vaf_df <- vaf_df[, c("Case", genesymbols_mut_driver, genesymbols_mut_othersmg)]

# make matrix for aa change-------------------------------------------------------------
aachange_mat <- get_somatic_mutation_detailed_matrix(pair_tab = genes2process, maf = maf_filtered_df)
aachange_df <- t(aachange_mat[,-1]) %>% as.data.frame()
## order
genesymbols_mut <- colnames(aachange_df)
aachange_df$Case <- rownames(aachange_df)
genesymbols_aachange_driver <- ccRCC_drivers[ccRCC_drivers %in% genesymbols_mut]
genesymbols_aachange_othersmg <- genesymbols_mut[!(genesymbols_mut %in% genesymbols_aachange_driver)]
aachange_df <- aachange_df[, c("Case", genesymbols_aachange_driver, genesymbols_aachange_othersmg)]

# make matrix for aa change- & vaf------------------------------------------------------------
mut_mat <- get_somatic_mutation_aachange_vaf_matrix(pair_tab = genes2process, maf = maf_filtered_df)
mut_df <- t(mut_mat[,-1]) %>% as.data.frame()
## order
genesymbols_mut <- colnames(mut_df)
mut_df$Case <- rownames(mut_df)
genesymbols_mut_driver <- ccRCC_drivers[ccRCC_drivers %in% genesymbols_mut]
genesymbols_mut_othersmg <- genesymbols_mut[!(genesymbols_mut %in% genesymbols_mut_driver)]
mut_df <- mut_df[, c("Case", genesymbols_mut_driver, genesymbols_mut_othersmg)]
mut_df <- mut_df %>%
  arrange(PBRM1, BAP1, SETD2, KDM5C)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "snRNA_ccRCC_Mutation_VAF_Table.", run_id, ".csv")
write.table(x = vaf_df, file = file2write, sep = ",", row.names = F, col.names = T, quote = F)
file2write <- paste0(dir_out, "snRNA_ccRCC_Mutation_AAchange_Table.", run_id, ".csv")
write.table(x = aachange_df, file = file2write, sep = ",", row.names = F, col.names = T, quote = F)
file2write <- paste0(dir_out, "snRNA_ccRCC_Mutation_AAchange_VAF_Table.", run_id, ".csv")
write.table(x = mut_df, file = file2write, sep = ",", row.names = F, col.names = T, quote = F)