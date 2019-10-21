# Yige Wu @ WashU 2019 Oct
## annotate the CPTAC3 ccRCC discovery set cases with mutation status
### mutation: VHL, PBRM1, BAP1, SETD2

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

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
snRNA_sample_mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Sample_Info/03_Post_request_snRNA_sample_selection/RCC_Specimen_Tracking - Case_Matrix.20191007.v1.tsv", data.table = F)

## 65 samples from 30 cases
case_annotation <- cptac_sample_map %>%
  filter(Type == "Tumor") %>%
  filter(Case.ID %in% snRNA_sample_mat$Case) %>%
  mutate(partID = Case.ID) %>%
  select(partID, Specimen.Label)

# input tumor purity ESTIMATE-RNA based result ---------------------------------------------------
# estimate_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S7.xlsx", sheet = "ESTIMATE scores")
# estimate_purity_rna <- as.numeric(as.vector(as.data.frame(estimate_tab[7,2:176])))
# estimate_purity_rna_df <- data.frame(ESTIMATE_TumorPurity_RNA = estimate_purity_rna, sampID = unlist(estimate_tab[3,2:176]))
# estimate_purity_rna_df$partID <- specimen2case_map[as.vector(estimate_purity_rna_df$sampID), "Case.ID"]
# estimate_purity_rna_df <- estimate_purity_rna_df %>%
#   dplyr::filter(!is.na(partID))
# 
# case_annotation <- merge(case_annotation, estimate_purity_rna_df %>%
#                              select(partID, ESTIMATE_TumorPurity_RNA), by = c("partID"), all.x = T)

# input mutation status ---------------------------------------------------
genes2add <- SMGs[["CCRCC"]]
maf_tab <- loadMaf()
mut_mat <- get_somatic_mutation_detailed_matrix(pair_tab = genes2add, maf = maf_tab)
mut2merge <- t(mut_mat) %>%
  as.data.frame()
mut2merge$partID <- rownames(mut2merge)
mut2merge %>% head()
case_annotation <- merge(case_annotation, mut2merge, by = c("partID"), all.x = T)

# input xcell result and add immnue group ------------------------------------------------------
xcell_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S7.xlsx", sheet = "xCell Signatures", skip = 2)
## immnue group 
immune_groups <- unlist(xcell_tab[1,2:ncol(xcell_tab)])
immune_group_df <- data.frame(immune_group = immune_groups, sampID = names(immune_groups))
case_annotation <- merge(case_annotation, immune_group_df, by.x = c("Specimen.Label"), by.y = c("sampID"), all.x = T)

# write out sample annotation table ---------------------------------------
case_annotation <- case_annotation %>%
  rename(Case = partID) %>%
  rename(Aliquot = Specimen.Label) %>%
  select(Case, Aliquot, VHL, PBRM1, BAP1, SETD2, MTOR, PTEN, TP53, KDM5C)
case_annotation_fn <- paste0(dir_out, "snRNA_ccRCC_Mutation_Table.", run_id, ".csv")
write.table(x = case_annotation, file = case_annotation_fn, sep = ",", row.names = F, col.names = T, quote = F)
