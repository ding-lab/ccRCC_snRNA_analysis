# Yige Wu @ WashU 2019 Aug
## annotate the CPTAC3 ccRCC discovery set cases with mutation status, subtype info and MET outlier status
### mutation: VHL, PBRM1, BAP1, SETD2
### subtypes: immune
### MET outlier status: protein outlier and phosphoprotein outlier
## add tumor Volume and normal Volume

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# set variables -----------------------------------------------------------
cancers2process <- c("CCRCC")

## variables for inputting outlier score table
outlier_sd <- 1.5

# input sample mapping file -----------------------------------------------
cptac_sample_map <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv", data.table = F)
cptac_sample_map <- data.frame(cptac_sample_map)
specimen2case_map <- cptac_sample_map %>%
  dplyr::filter(Type == "Tumor") %>%
  select(Case.ID, Specimen.Label)
rownames(specimen2case_map) <- specimen2case_map$Specimen.Label

cptac_shiped2wustl_sample_map <- read_excel("Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/sample_info/CPTAC3_ccRCC_request/CPTAC3_inventory/JHU Distribution Manifest_ccRCC material to WUSTL_6.17.2019.xlsx")
cptac_shiped2wustl_sample_map %>%
  select(`Subject ID`) %>%
  unique() %>%
  nrow()
## 65 samples from 30 cases
case_annotation <- cptac_sample_map %>%
  dplyr::filter(Type == "Tumor") %>%
  dplyr::filter(Case.ID %in% cptac_shiped2wustl_sample_map$`Subject ID`) %>%
  mutate(partID = Case.ID) %>%
  select(partID, Specimen.Label)

# input tumor purity ESTIMATE-RNA based result ---------------------------------------------------
estimate_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S7.xlsx", sheet = "ESTIMATE scores")
estimate_purity_rna <- as.numeric(as.vector(as.data.frame(estimate_tab[7,2:176])))
estimate_purity_rna_df <- data.frame(ESTIMATE_TumorPurity_RNA = estimate_purity_rna, sampID = unlist(estimate_tab[3,2:176]))
estimate_purity_rna_df$partID <- specimen2case_map[as.vector(estimate_purity_rna_df$sampID), "Case.ID"]
estimate_purity_rna_df <- estimate_purity_rna_df %>%
  dplyr::filter(!is.na(partID))

case_annotation <- merge(case_annotation, estimate_purity_rna_df %>%
                             select(partID, ESTIMATE_TumorPurity_RNA), by = c("partID"), all.x = T)

# input mutation status ---------------------------------------------------
genes2add <- c("VHL", "PBRM1", "BAP1", "SETD2")
maf_tab <- loadMaf()
mut_mat <- generate_somatic_mutation_matrix(pair_tab = genes2add, maf = maf_tab)
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

# add MET protein outlier data ------------------------------------------------------
partIDs <- case_annotation$partID
file2input <- paste0("./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/escore_tab_", "CCRCC", "_", "kinase", "_reg_nonNA", "20", ".txt")
escore_tab <- fread(input = file2input, data.table = F)
escore_tab_scaled <- escore_tab[, partIDs]
escore_tab_scaled <- scale_by_row(escore_tab_scaled)
escore_tab_outlier <- cbind(escore_tab[, c("pair", "SELF", "GENE")], (escore_tab_scaled > outlier_sd))

pro_outlier_tmp <- escore_tab_outlier[escore_tab_outlier$SELF == "cis" & escore_tab_outlier$GENE == "MET",]
pro_outlier_tmp <- pro_outlier_tmp[1,]
if (!is.null(pro_outlier_tmp)) {
  if (nrow(pro_outlier_tmp) > 0) {
    pro_outlier.m <- melt(pro_outlier_tmp, 
                          id.vars = colnames(pro_outlier_tmp)[!(colnames(pro_outlier_tmp) %in% partIDs)])
    pro_outlier2merge <- pro_outlier.m %>%
      mutate(partID = variable) %>%
      mutate(MET_PRO = value) %>%
      select(partID, MET_PRO)
    case_annotation <- merge(case_annotation, pro_outlier2merge, by = c("partID"), all.x = T)
  }
}

# add MET Y1234 phosphorylation outlier -----------------------------------------------
file2input <- paste0("./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/sscore_tab_", "CCRCC", "_", "kinase", "_reg_nonNA", "20", ".txt")
sscore_tab <- fread(input = file2input, data.table = F)
sscore_tab_scaled <- sscore_tab[, partIDs]
sscore_tab_scaled <- scale_by_row(sscore_tab_scaled)
sscore_tab_outlier <- cbind(sscore_tab[, c("SUB_GENE", "SUB_MOD_RSD")], (sscore_tab_scaled > outlier_sd))


phog_outlier_tmp <- sscore_tab_outlier[sscore_tab_outlier$SUB_MOD_RSD == "Y1234" & sscore_tab_outlier$SUB_GENE == "MET",]
phog_outlier_tmp <- phog_outlier_tmp[1,]
if (!is.null(phog_outlier_tmp)) {
  if (nrow(phog_outlier_tmp) > 0) {
    phog_outlier.m <- melt(phog_outlier_tmp, 
                           id.vars = colnames(phog_outlier_tmp)[!(colnames(phog_outlier_tmp) %in% partIDs)])
    phog_outlier2merge <- phog_outlier.m %>%
      mutate(partID = variable) %>%
      mutate(MET_Y1234_PHO = value) %>%
      select(partID, MET_Y1234_PHO)
    case_annotation <- merge(case_annotation, phog_outlier2merge, by = c("partID"), all.x = T)
  }
}


# add the number of tumor segment and the Volume of tumor -----------------------------------------
tumor_Volume_segments_tab <- read_excel(path = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/sample_info/CPTAC3_ccRCC_request/CPTAC3_inventory/JHU Distribution Manifest_ccRCC material to WUSTL_6.17.2019.xlsx")
tumor_Volume_frozen_powder_segments_tab <- tumor_Volume_segments_tab %>%
  dplyr::filter(Preparation == "Frozen") %>%
  dplyr::filter(`Material Type` == "pulverized tissue")
case_annotation$Num_Frozen_pulverized_Tumor_segments <- sapply(X = case_annotation$partID, FUN = function(partID, tumor_Volume_segments_tab) nrow(tumor_Volume_segments_tab[tumor_Volume_segments_tab$`Subject ID` == partID,]), tumor_Volume_segments_tab = tumor_Volume_frozen_powder_segments_tab)
case_annotation$Volume_Frozen_pulverized_Tumor_mg <-  sapply(X = case_annotation$partID, FUN = function(partID, tumor_Volume_segments_tab) paste0(tumor_Volume_segments_tab$Volume[tumor_Volume_segments_tab$`Subject ID` == partID], collapse = "|"), tumor_Volume_segments_tab = tumor_Volume_frozen_powder_segments_tab)

# add the Volume of NAT ---------------------------------------------------
NAT_Volume_segments_tab <- read_excel(path = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/sample_info/CPTAC3_ccRCC_request/selected_CCRCC_samples_v7.xlsx")
NAT_Volume_segments_tab <- data.frame(NAT_Volume_segments_tab)
NAT_Volume2merge <- NAT_Volume_segments_tab %>%
  dplyr::select(partID, Volume_NAT)
case_annotation <- merge(case_annotation, NAT_Volume2merge, all.x = T)
# write out sample annotation table ---------------------------------------
case_annotation_fn <- paste0(makeOutDir(), "ccRCC_snRNA_case_annotation_v2.csv")
write.table(x = case_annotation, file = case_annotation_fn, sep = ",", row.names = F, col.names = T, quote = F)
