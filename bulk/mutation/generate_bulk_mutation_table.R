# Yige Wu @ WashU 2019 Oct
## annotate the CPTAC3 ccRCC discovery set cases with mutation status
### mutation: VHL, PBRM1, BAP1, SETD2

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

# functions ---------------------------------------------------------------
get_somatic_mutation_vaf_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  maf$vaf <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    VAF <- paste0(unique(x), collapse = ",")
    return(VAF)
  }, value.var = "vaf", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

get_somatic_mutation_detailed_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  }, value.var = "HGVSp_Short", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}


# input dependencies ------------------------------------------------------
## input maf file
maf_df <- loadMaf()
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# set parameters ----------------------------------------------------------
## set genes to look for mutations
genes2process <- ccRCC_SMGs
## set sample case ids to process
caseids2process <- unique(idmetadata_df$Case[idmetadata_df$snRNA_available])
## filter maf file
maf_df <- 

# make matrix for vaf-------------------------------------------------------------
mut_mat <- get_somatic_mutation_vaf_matrix(pair_tab = genes2process, maf = maf_df)
mut_t_df <- t(mut_mat[,-1]) %>% as.data.frame()
mut2merge <- t(mut_mat) %>%
  as.data.frame()
mut2merge$partID <- rownames(mut2merge)
mut2merge %>% head()
case_annotation <- merge(case_annotation, mut2merge, by = c("partID"), all.x = T)

case_annotation <- case_annotation %>%
  rename(Case = partID) %>%
  rename(Aliquot = Specimen.Label) %>%
  select(Case, Aliquot, VHL, PBRM1, BAP1, SETD2, MTOR, PTEN, TP53, KDM5C)


# write output ------------------------------------------------------------
case_annotation_fn <- paste0(dir_out, "snRNA_ccRCC_Mutation_VAF_Table.", run_id, ".csv")
write.table(x = case_annotation, file = case_annotation_fn, sep = ",", row.names = F, col.names = T, quote = F)
