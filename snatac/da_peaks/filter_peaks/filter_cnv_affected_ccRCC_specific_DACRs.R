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
dacrs_up_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/da_up_peaks_Tumor_vs_PT.annotated.Included_all_nFC.20210603.tsv")
dacrs_down_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/da_down_peaks_Tumor_vs_PT.annotated.Included_all_nFC.20210603.tsv")

# preprocess ----------------------------
## get the aliquot ids to process
cases_process <- unique(metadata_df$Case[metadata_df$snATAC_available & metadata_df$snATAC_used])
## count the number of up/down fold changes including the insignificant fold changes
colnames_allfcs <- colnames(dacrs_up_df)[grepl(x = colnames(dacrs_up_df), pattern = "All_avg_lnFC")]
dacrs_up_df$Count_oppositeFC <- rowSums(dacrs_up_df[,colnames_allfcs] < 0, na.rm = T)
colnames_allfcs <- colnames(dacrs_down_df)[grepl(x = colnames(dacrs_down_df), pattern = "All_avg_lnFC")]
dacrs_down_df$Count_oppositeFC <- rowSums(dacrs_down_df[,colnames_allfcs] > 0, na.rm = T)
## get all the types of DACR gene annotation type
dacrs_df <- rbind(dacrs_up_df %>%
                    mutate(DAP_type = str_split_fixed(string = Type, pattern = " \\(", n = 2)[,1]) %>%
                    mutate(DAP_direction = "Up") %>%
                    rename(Count_sig = Count_up) %>%
                    select(peak, Gene, DAP_type, DAP_direction, Count_sig, Count_oppositeFC),
                  dacrs_down_df %>%
                    mutate(DAP_type = str_split_fixed(string = Type, pattern = " \\(", n = 2)[,1]) %>%
                    mutate(DAP_direction = "Down") %>%
                    rename(Count_sig = Count_down) %>%
                    select(peak, Gene, DAP_type, DAP_direction, Count_sig, Count_oppositeFC))
table(dacrs_df$DAP_type)
## get genes to process
genes_process <- unique(dacrs_df$Gene)

# preprocess mean CNV values per gene--------------------------------------------------------------
## filter the CNVs
cna_filtered_df <- cna_df[cna_df$gene_name %in% genes_process,]
## preprocess the CNV data frame
colnames_old <- colnames(cna_filtered_df)
colnames_new <- str_split_fixed(string = colnames_old, pattern = "\\.", n = 4)[,1]
colnames(cna_filtered_df) <- colnames_new
rownames(cna_filtered_df) <- cna_filtered_df$gene_name
## filter columns
cna_filtered_df <- cna_filtered_df[, c("gene_id", "gene_name", cases_process)]
## calculate mean CNV value
cna_filtered_df$cnv_mean <- rowMeans(cna_filtered_df[, cases_process])
summary(cna_filtered_df$cnv_mean)
## map back to peaks
dacrs_df$gene_cnv_mean <- mapvalues(x = dacrs_df$Gene, from = cna_filtered_df$gene_name, to = as.vector(cna_filtered_df$cnv_mean))
dacrs_df$gene_cnv_mean[dacrs_df$gene_cnv_mean == dacrs_df$Gene] <- NA
dacrs_df$gene_cnv_mean <- as.numeric(as.vector(dacrs_df$gene_cnv_mean))
## decide potential CNv effect
dacrs_df <- dacrs_df %>%
  mutate(potential_cnv_effect = ifelse(gene_cnv_mean >= 0.1, "Up",
                                       ifelse(gene_cnv_mean <= -0.1, "Down", "Not sure"))) %>%
  mutate(chr_peak = str_split_fixed(string = peak, pattern = "\\-", n = 3)[,1])
length(which(is.na(dacrs_df$gene_cnv_mean))) # [1] 5230
length(which(is.na(dacrs_df$gene_cnv_mean)))/nrow(dacrs_df) # [1] 0.189816
## filter
dacrs_filtered_df <- dacrs_df %>%
  filter(Count_sig >= 12) %>%
  filter(Count_oppositeFC == 0) %>%
  filter(!is.na(potential_cnv_effect) & potential_cnv_effect != DAP_direction)

# test effect of the removal ----------------------------------------------
dacrs_df %>%
  filter(Count_sig >= 12) %>%
  filter(Count_oppositeFC == 0) %>%
  select(DAP_direction, chr_peak) %>%
  table()
dacrs_filtered_df %>%
  filter(Count_sig >= 12) %>%
  filter(Count_oppositeFC == 0) %>%
  select(DAP_direction, chr_peak) %>%
  table()

dacrs_df %>%
  filter(Count_sig >= 12) %>%
  filter(Count_oppositeFC == 0) %>%
  nrow()

dacrs_filtered_df %>%
  filter(Count_sig >= 12) %>%
  nrow()

# write table -------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_specific.DACRs.PotentialCNVEffectAnnotated.", run_id, ".tsv")
write.table(x = dacrs_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ccRCC_specific.DACRs.PotentialCNVEffectFiltered.", run_id, ".tsv")
write.table(x = dacrs_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
