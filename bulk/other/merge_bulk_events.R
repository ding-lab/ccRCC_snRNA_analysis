# Yige Wu @WashU March 2020
## running on local
## for making a table including mutation (SMG), CNV (reported by TCGA, bulk), VHL methylation and SV info detected by bulk omics by sample
## By sample

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input bulk mutation result
maf_df <- loadMaf()
## input bulk CNV profile
arm_cnv_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/bulk/copy_number/write_sample_bicseq_cnv_profile/20200316.v1/Bulk_WGS_Chr_CNV_Profile.20200316.v1.tsv", data.table = F)
## input methylatin value
methylation_by_gene_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/bulk/methylation/average_methylation_by_gene_by_aliquot/20200317.v1/CCRCC.TumorBeta.by_gene.20200317.v1.tsv", data.table = F)
## input chr3 translocation data
chr3_translocation_df <- openxlsx::read.xlsx("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/CPTAC-3-ccRCC_paper_STable S2.xlsx", sheet = "Tab1.Chr3_translocation_all", fillMergedCells = T)
## input case ids to process
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_cluster/write_path_to_seurat_objects_on_box/20200219.v1/Seurat_Object_Paths.20200219.v1.tsv", data.table = F)

# format all the dataset ready for merging --------------------------------
## set case ids to process
case_ids <- unique(srat_paths$Case)

## format mutation data
maf_df <- maf_df %>%
  mutate(case_id = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", n = 2)[,1]) %>%
  filter(case_id %in% case_ids)
var_class_df <- generate_somatic_mutation_matrix(pair_tab = ccRCC_SMGs, maf = maf_df)
var_class_by_case_mat <- t(var_class_df[,-1])
var_class_by_case_mat[var_class_by_case_mat == "Silent"] <- ""
var_class_by_case_df <- as.data.frame(var_class_by_case_mat)
### change column names
colnames(var_class_by_case_df) <- paste0("Mut.", colnames(var_class_by_case_df))
var_class_by_case_df$Case <- rownames(var_class_by_case_df)

## format arm level cnv
arm_cnv_by_case_df <- arm_cnv_df %>%
  dplyr::select("Case", 'CN.3p', 'CN.5q', "CN.14q") %>%
  dplyr::filter(Case %in% case_ids)

## format translocation data
### get the unique chr3 translocation for each case
uniq_chr3_translocation_long_df <- chr3_translocation_df %>%
  select(CASE_ID, TYPE) %>%
  unique() %>%
  as.data.frame() %>%
  rename(Case = CASE_ID) %>%
  rename(Chr3_Translocation_Type = TYPE) %>%
  filter(Case %in% case_ids) %>%
  arrange(Chr3_Translocation_Type)
uniq_chr3_translocation_long_df1 <- uniq_chr3_translocation_long_df %>%
  filter(!duplicated(Case))
uniq_chr3_translocation_long_df2 <- uniq_chr3_translocation_long_df %>%
  filter(duplicated(Case))
### merge the first translocation with the second
chr3_translocation_wide_df <- merge(uniq_chr3_translocation_long_df1, uniq_chr3_translocation_long_df2, by = c("Case"), suffixes = c("1", "2"), all = T)

## format methylation data
### filter down to VHL only
methylation_vhl_wide_df <- methylation_by_gene_df %>%
  filter(gene_symbol == "VHL")
methylation_vhl_wide_mat <- methylation_vhl_wide_df[,-1]
methylation_vhl_long_df <- t(methylation_vhl_wide_mat)
colnames(methylation_vhl_long_df) <- "Methyl.VHL"
methylation_vhl_long_df <- as.data.frame(methylation_vhl_long_df)
methylation_vhl_long_df$Case <- rownames(methylation_vhl_long_df)
# merge data --------------------------------------------------------------
merged_bulk_events_df <- data.frame(Case = case_ids)
merged_bulk_events_df <- merge(merged_bulk_events_df, var_class_by_case_df, by = c("Case"), all.x = T)
merged_bulk_events_df <- merge(merged_bulk_events_df, arm_cnv_by_case_df, by = c("Case"), all.x = T)
merged_bulk_events_df <- merge(merged_bulk_events_df, chr3_translocation_wide_df, by = c("Case"), all.x = T)
merged_bulk_events_df <- merge(merged_bulk_events_df, methylation_vhl_long_df, by = c("Case"), all.x = T)

# write table -------------------------------------------------------------
file2write <- paste0(dir_out, "merged_bulk_events.", run_id, ".tsv")
write.table(x = merged_bulk_events_df, file = file2write, sep = "\t", quote = F, row.names = F)

