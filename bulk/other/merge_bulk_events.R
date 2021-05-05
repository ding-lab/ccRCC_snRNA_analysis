# Yige Wu @WashU March 2020
## running on local
## for making a table including mutation (SMG), CNV (reported by TCGA, bulk), VHL methylation and SV info detected by bulk omics by sample
## By sample

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
## input bulk mutation result
maf_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/merge_discovery_ith_mutations_for_snRNA_samples/20210504.v1/snRNA_samples.somatic_mutation.maf.v1.0.tsv")
## input germline mutation result
vhl_germline_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Germline_Variants/Analysis_Results/2019-09-10/VHL_var.txt")
bap1_germline_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Germline_Variants/Analysis_Results/2019-09-10/BAP1_var.txt")
## input bulk CNV profile
arm_cnv_df <- fread(input = "./Resources/Analysis_Results/bulk/copy_number/write_sample_bicseq_cnv_profile/20210503.v1/Bulk_WGS_Chr_CNV_Profile.20210503.v1.tsv", data.table = F)
## input methylatin value
methylation_by_gene_df <- fread(input = "./Resources/Analysis_Results/bulk/methylation/average_methylation_by_gene_by_aliquot/20200317.v1/CCRCC.TumorBeta.by_gene.20200317.v1.tsv", data.table = F)
## input chr3 translocation data
chr3_translocation_df <- openxlsx::read.xlsx("./Resources/Bulk_Processed_Data/CPTAC3-ccRCC-SupplementaryTables_Final/CPTAC-3-ccRCC_paper_STable S2.xlsx", sheet = "Tab1.Chr3_translocation_all", fillMergedCells = T)
## input id meta data
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv", data.table = F)
## input 10xmapping results
sn_mut_df <- fread(input = "./Resources/Analysis_Results/mutation/unite_10xmapping/20200303.v1/10XMapping.20200303.v1.tsv", data.table = F)


# preprocess --------------------------------------------------------------
## set case ids to process
case_ids <- unique(id_metadata_df$Case[id_metadata_df$snRNA_available])
sample_ids <- unique(id_metadata_df$Sample[id_metadata_df$snRNA_available])
genes_mut <-genes_pathogenicpathways_in_ccRCC

# format somatic mutation data for the discovery data ----------------------------------------------------
get_mutation_class_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat & maf$Variant_Classification != "Silent",]
  nrow(maf)
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sample_id, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  }, value.var = "Variant_Classification", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}
var_class_df <- get_mutation_class_matrix(pair_tab = genes_mut, maf = maf_df)
var_class_by_case_mat <- t(var_class_df[,-1])
var_class_by_case_mat[var_class_by_case_mat == "Silent"] <- ""
var_class_by_case_mat[var_class_by_case_mat == ""] <- "None"
var_class_by_case_df <- as.data.frame(var_class_by_case_mat)
### change column names
colnames(var_class_by_case_df) <- paste0("Mut.", colnames(var_class_by_case_df))
var_class_by_case_df$Sample <- rownames(var_class_by_case_df)
var_class_by_case_df$Case <- mapvalues(x = var_class_by_case_df$Sample, from = id_metadata_df$Sample, to = as.vector(id_metadata_df$Case))
## add new aliquot id
var_class_by_case_df$Aliquot_snRNA_WU <- mapvalues(x = var_class_by_case_df$Sample, from = id_metadata_df$Sample, to = as.vector(id_metadata_df$Aliquot.snRNA.WU))

# add germline mutation data ----------------------------------------------
var_class_by_case_df[, "Mut.VHL.Germline"] <- "None"
var_class_by_case_df[var_class_by_case_df$Aliquot_snRNA_WU == "C3L-00917-T1", "Mut.VHL.Germline"] <- "Nonsense_Mutation"
var_class_by_case_df[, "Mut.BAP1.Germline"] <- "None"
var_class_by_case_df[var_class_by_case_df$Aliquot_snRNA_WU == "C3N-00177-T1", "Mut.BAP1.Germline"] <- "Frame_Shift_Del"

# add 10xmapping data -----------------------------------------------------
sn_mut_df$Aliquot_snRNA_WU <- mapvalues(x = sn_mut_df$aliquot, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Aliquot.snRNA.WU))
sn_mut_df$Case <- mapvalues(x = sn_mut_df$aliquot, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))
sn_mut_df <- sn_mut_df %>%
  filter(allele_type == "Var") %>%
  filter(gene_symbol %in% var_class_df$Hugo_Symbol) %>%
  # filter(!(grepl(x = Aliquot_snRNA_WU, pattern = "T1"))) %>%
  filter(!(Aliquot_snRNA_WU %in% var_class_by_case_df$Aliquot_snRNA_WU))
## loop to add the mutation info for the additional segment
var_class_additional_df <- NULL
for (aliquot_wu_tmp in unique(sn_mut_df$Aliquot_snRNA_WU)) {
  case_tmp <- sn_mut_df$Case[sn_mut_df$Aliquot_snRNA_WU == aliquot_wu_tmp]
  var_tmp <- var_class_by_case_df %>%
    filter(Case == case_tmp)
  ## change based on the discovery mutation info
  var_tmp_add <- var_tmp[1,]
  var_tmp_add$Aliquot_snRNA_WU <- aliquot_wu_tmp
  var_tmp_add$Sample <- mapvalues(x = aliquot_wu_tmp, from = id_metadata_df$Aliquot.snRNA.WU, to = as.vector(id_metadata_df$Sample))
  
  sn_mut_tmp <- sn_mut_df %>%
    filter(Aliquot_snRNA_WU == aliquot_wu_tmp)
  genes_mut_detected <- unique(sn_mut_tmp$gene_symbol)
  ## see if any mutation not detected
  genes_mut <- colnames(var_tmp)
  genes_mut <- genes_mut[grepl(x = genes_mut, pattern = "Mut")]
  genes_mut <- gsub(x = genes_mut, pattern = "Mut.", replacement = "")
  genes_mut_notdetected <- genes_mut[!(genes_mut %in% genes_mut_detected)]
  if (length(genes_mut_notdetected) > 0) {
    ### any undetected mutation, change the mutation status to NA: data unavailable/unknown
    var_tmp_add[, paste0("Mut.",  genes_mut_notdetected)] <- NA
  }
  var_tmp_add[var_tmp_add == "None"] <- NA
  var_class_additional_df <- rbind(var_class_additional_df, var_tmp_add)
}
var_class_df <- rbind(var_class_by_case_df, var_class_additional_df)

# format arm level cnv ----------------------------------------------------
arm_cnv_by_case_df <- arm_cnv_df %>%
  dplyr::select("Case", 'CN.3p', 'CN.5q', "CN.14q") %>%
  rename(CN.bulk.3p = CN.3p) %>%
  rename(CN.bulk.5q = CN.5q) %>%
  rename(CN.bulk.14q = CN.14q) %>%
  dplyr::filter(Case %in% case_ids)
arm_cnv_by_case_df[arm_cnv_by_case_df == ""] <- "neutral"
arm_cnv_by_case_df$Aliquot_snRNA_WU <- paste0(arm_cnv_by_case_df$Case, "-T1")

# format translocation data -----------------------------------------------
### get the unique chr3 translocation for each case
uniq_chr3_translocation_long_df <- chr3_translocation_df %>%
  select(CASE_ID, TYPE) %>%
  unique() %>%
  as.data.frame() %>%
  rename(Case = CASE_ID) %>%
  rename(Translocation = TYPE) %>%
  filter(Case %in% case_ids) %>%
  arrange(Translocation)
chr35_translocation_long_df <- uniq_chr3_translocation_long_df %>%
  filter(Translocation == "chr3-5")
chr32_translocation_long_df <- uniq_chr3_translocation_long_df %>%
  filter(Translocation == "chr3-2")
chr3_translocation_other_long_df <- uniq_chr3_translocation_long_df %>%
  filter(!(Translocation %in% c("chr3-2", "chr3-5"))) %>%
  rename(Translocation.t3_other = Translocation)
### merge the first translocation with the second
chr3_translocation_wide_df <- merge(chr35_translocation_long_df, chr32_translocation_long_df, by = c("Case"), suffixes = c(".t35", ".t32"), all = T)
chr3_translocation_wide_df <- merge(chr3_translocation_wide_df, chr3_translocation_other_long_df, by = c("Case"), all = T)
chr3_translocation_wide_df$Aliquot_snRNA_WU <- paste0(chr3_translocation_wide_df$Case, "-T1")

# format methylation data -------------------------------------------------
### filter down to VHL only
methylation_vhl_wide_df <- methylation_by_gene_df %>%
  filter(gene_symbol == "VHL")
methylation_vhl_wide_mat <- methylation_vhl_wide_df[,-1]
methylation_vhl_long_df <- t(methylation_vhl_wide_mat)
colnames(methylation_vhl_long_df) <- "Methyl.VHL"
methylation_vhl_long_df <- as.data.frame(methylation_vhl_long_df)
methylation_vhl_long_df$Case <- rownames(methylation_vhl_long_df)
methylation_vhl_long_df$Aliquot_snRNA_WU <- paste0(methylation_vhl_long_df$Case, "-T1")

# merge data --------------------------------------------------------------
merged_bulk_events_df <- data.frame(Case = case_ids)
merged_bulk_events_df <- merge(merged_bulk_events_df, var_class_df, by = c("Case"), all = T)
merged_bulk_events_df <- merge(merged_bulk_events_df, arm_cnv_by_case_df, by = c("Case", "Aliquot_snRNA_WU"), all.x = T)
merged_bulk_events_df <- merge(merged_bulk_events_df, chr3_translocation_wide_df, by = c("Case", "Aliquot_snRNA_WU"), all.x = T)
## fill in empty cells without events
merged_bulk_events_T1_df <- merged_bulk_events_df[grepl(x = merged_bulk_events_df$Aliquot_snRNA_WU, pattern = "T1"),]
merged_bulk_events_T23_df <- merged_bulk_events_df[!grepl(x = merged_bulk_events_df$Aliquot_snRNA_WU, pattern = "T1"),]
merged_bulk_events_T1_df[is.na(merged_bulk_events_T1_df)] <- "None"
merged_bulk_events_df <- rbind(merged_bulk_events_T1_df, merged_bulk_events_T23_df)
## fill in the methylation data
merged_bulk_events_df <- merge(merged_bulk_events_df, methylation_vhl_long_df, by = c("Case", "Aliquot_snRNA_WU"), all.x = T)

# write table -------------------------------------------------------------
file2write <- paste0(dir_out, "merged_bulk_events.", run_id, ".tsv")
write.table(x = merged_bulk_events_df, file = file2write, sep = "\t", quote = F, row.names = F)

