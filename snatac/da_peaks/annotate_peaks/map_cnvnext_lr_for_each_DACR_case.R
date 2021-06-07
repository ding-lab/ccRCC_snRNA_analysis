# Yige Wu @ WashU 2021 Jun
## annotate sample copy number profile (3p, 5q, 14q)
## CNV: https://wustl.box.com/s/vlde6w791k81q1aibvgs0a487zxleij6
# Purity and ploidy: https://wustl.box.com/s/jel5krgvnvlq5z32vdg4itdcq6znuqzn
# From UMich

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
## input CNA matrix
cna_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Combined/Absolute_cnv/c3-ccrcc-combined-cnvex-lr_v1.0.csv")
## input snRNA sample set
metadata_df <- fread("./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv", data.table = F)
## input DACRs
dacrs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/filter_peaks/filter_cnv_affected_ccRCC_specific_DACRs/20210603.v1/ccRCC_specific.DACRs.PotentialCNVEffectAnnotated.20210603.v1.tsv")

# preprocess ----------------------------
## get genes to process
dacrs_filtered_df <- dacrs_df %>%
  filter(Count_sig >= 12) %>%
  filter(Count_oppositeFC == 0)
genes_process <- unique(dacrs_filtered_df$Gene)
## get cases to process
cases_process <- metadata_df$Case[metadata_df$snATAC_available & metadata_df$snATAC_used]
cases_process <- unique(cases_process)

# preprocess mean CNV values per gene--------------------------------------------------------------
## filter the CNVs
cna_filtered_df <- cna_df[cna_df$gene_name %in% genes_process,]
## preprocess the CNV data frame
colnames_old <- colnames(cna_filtered_df)
colnames_new <- str_split_fixed(string = colnames_old, pattern = "\\.", n = 4)[,1]
colnames(cna_filtered_df) <- colnames_new
cna_filtered_df <- cna_filtered_df[!duplicated(cna_filtered_df$gene_name),]
## add peak
cna_bypeak_df <- merge(x = cna_filtered_df, 
                       y = dacrs_filtered_df %>%
                         select(peak, Gene),
                       by.x = c("gene_name"), by.y = c("Gene"), all.x = T)

rownames(cna_bypeak_df) <- cna_bypeak_df$peak
cna_t_mat <- t(as.matrix(cna_bypeak_df[,cases_process]))
cna_t_mat[1:5, 1:5]

# write table -------------------------------------------------------------
file2write <- paste0(dir_out, "Case2Peak.CNV.", run_id, ".RDS")
saveRDS(object = cna_t_mat, file = file2write, compress = T)
