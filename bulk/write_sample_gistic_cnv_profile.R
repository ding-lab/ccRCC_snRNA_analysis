# Yige Wu @ WashU 2020 Feb
## annotate sample copy number profile (3p, 5q, 14q)

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
## input the list of cases with snRNA data
srat_paths <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_cluster/write_path_to_seurat_objects_on_box/20200219.v1/Seurat_Object_Paths.20200219.v1.tsv", data.table = F)
## input 3p CNV profile
cnv_3p_df <- read_excel(path = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/CPTAC-3-ccRCC_paper_STable S2.xlsx", sheet = "Tab5.3p_events")
## input 5q CNV profile
cnv_5q_df <- read_excel(path = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/CPTAC-3-ccRCC_paper_STable S2.xlsx", sheet = "Tab6.5q_events")
## input 14q CNV profile
cnv_14q_df <- read_excel(path = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/CPTAC-3-ccRCC_paper_STable S2.xlsx", sheet = "Tab8.14q_events")

# merge -------------------------------------------------------------------
chr_cnv_state_df <- data.frame(Case = unique(srat_paths$Case))
## merge with 3p CNV
colnames(cnv_3p_df) <- c("Case", "Aliquot.DNA", "GISTIC2.3p", "Pathologist_Review.3p")
chr_cnv_state_df <- merge(chr_cnv_state_df, cnv_3p_df, by = c("Case"), all.x = T)
## merge with 5q CNV
colnames(cnv_5q_df) <- c("Case", "Aliquot.DNA", "GISTIC2.5q", "Pathologist_Review.5q")
chr_cnv_state_df <- merge(chr_cnv_state_df, cnv_5q_df, by = c("Case", "Aliquot.DNA"),all.x = T)
## merge with 14q CNV
colnames(cnv_14q_df) <- c("Case", "Aliquot.DNA", "GISTIC2.14q", "Pathologist_Review.14q", "notes.14q")
chr_cnv_state_df <- merge(chr_cnv_state_df, cnv_14q_df, by = c("Case", "Aliquot.DNA"),all.x = T)

# write table -------------------------------------------------------------
write.table(x = chr_cnv_state_df, file = paste0(dir_out, "Bulk_WGS_GISTIC2_Chr_CNV_Profile", ".", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

