# Yige Wu @ WashU 2019 Oct
## annotate the CPTAC3 ccRCC discovery set cases with mutation status
### mutation: VHL, PBRM1, BAP1, SETD2

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


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


# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# set variables -----------------------------------------------------------
cancers2process <- c("CCRCC")

# input sample mapping file -----------------------------------------------
cptac_sample_map <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv", data.table = F)
cptac_sample_map <- data.frame(cptac_sample_map)

specimen2case_map <- cptac_sample_map %>%
  dplyr::filter(Type == "Tumor") %>%
  select(Case.ID, Specimen.Label)
rownames(specimen2case_map) <- specimen2case_map$Specimen.Label


# input the snRNA sample matrix -------------------------------------------
snRNA_sample_mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Sample_Info/03_Post_request_snRNA_sample_selection/RCC_Specimen_Tracking - Case_Matrix.20191007.v1.tsv", data.table = F)

# write out mutation short amino acid change matrix ---------------------------------------
case_annotation <- cptac_sample_map %>%
  filter(Type == "Tumor") %>%
  filter(Case.ID %in% snRNA_sample_mat$Case) %>%
  mutate(partID = Case.ID) %>%
  select(partID, Specimen.Label)

genes2add <- SMGs[["CCRCC"]]
maf_tab <- loadMaf()
mut_mat <- get_somatic_mutation_detailed_matrix(pair_tab = genes2add, maf = maf_tab)
mut2merge <- t(mut_mat) %>%
  as.data.frame()
mut2merge$partID <- rownames(mut2merge)
mut2merge %>% head()
case_annotation <- merge(case_annotation, mut2merge, by = c("partID"), all.x = T)

case_annotation <- case_annotation %>%
  rename(Case = partID) %>%
  rename(Aliquot = Specimen.Label) %>%
  select(Case, Aliquot, VHL, PBRM1, BAP1, SETD2, MTOR, PTEN, TP53, KDM5C)
case_annotation_fn <- paste0(dir_out, "snRNA_ccRCC_Mutation_Table.", run_id, ".csv")
write.table(x = case_annotation, file = case_annotation_fn, sep = ",", row.names = F, col.names = T, quote = F)


# write out mutation VAF matrix -------------------------------------------
case_annotation <- cptac_sample_map %>%
  filter(Type == "Tumor") %>%
  filter(Case.ID %in% snRNA_sample_mat$Case) %>%
  mutate(partID = Case.ID) %>%
  select(partID, Specimen.Label)

genes2add <- SMGs[["CCRCC"]]
maf_tab <- loadMaf()
mut_mat <- get_somatic_mutation_vaf_matrix(pair_tab = genes2add, maf = maf_tab)
mut2merge <- t(mut_mat) %>%
  as.data.frame()
mut2merge$partID <- rownames(mut2merge)
mut2merge %>% head()
case_annotation <- merge(case_annotation, mut2merge, by = c("partID"), all.x = T)

case_annotation <- case_annotation %>%
  rename(Case = partID) %>%
  rename(Aliquot = Specimen.Label) %>%
  select(Case, Aliquot, VHL, PBRM1, BAP1, SETD2, MTOR, PTEN, TP53, KDM5C)
case_annotation_fn <- paste0(dir_out, "snRNA_ccRCC_Mutation_VAF_Table.", run_id, ".csv")
write.table(x = case_annotation, file = case_annotation_fn, sep = ",", row.names = F, col.names = T, quote = F)
