# Yige Wu @WashU May 2021

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
coaccess_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Peak_Annotation/da_up_peaks_Tumor_vs_PT.annotated.CICERO.20210513.tsv")

# filter and extract ------------------------------------------------------
nrow(coaccess_df)
coaccess_filtered_df <- coaccess_df %>%
  filter(Count_up >= 5) %>%
  filter(!is.na(cicero_coaccess) & cicero_coaccess >= 0.25) %>%
  mutate(peak_coacess = str_split_fixed(string = cicero_co_accessible_peaks, pattern = "_", n = 2)[,2])
nrow(coaccess_filtered_df)
## first scenerio: DAP is promoter region for gene X, then the coaccessible peak is potential enhancer for gene X
coaccess_enhancer_df <- coaccess_filtered_df %>%
  filter(Type == "Promoter") %>%
  select(peak_coacess, Gene) %>%
  rename(Peak = peak_coacess) %>%
  mutate(Peak_Type = "Enhancer") %>%
  mutate(Is.CAP = T) %>%
  mutate(Is.DAP = F)
coaccess_enhancer_df$Is.DAP[coaccess_enhancer_df$Peak %in% coaccess_filtered_df$peak] <- T
## second scenerio: the coaccessible peak is promoter region for gene Y, then the DAP is potential enhancer for gene Y
coaccess_promoter_df <- coaccess_filtered_df %>%
  filter(cicero_Type == "Promoter") %>%
  select(peak_coacess, cicero_Gene, cicero_Type) %>%
  rename(Peak = peak_coacess) %>%
  rename(Gene = cicero_Gene) %>%
  rename(Peak_Type = cicero_Type) %>%
  mutate(Is.CAP = T) %>%
  mutate(Is.DAP = F)
dap_enhancer_df <- coaccess_filtered_df %>%
  filter(cicero_Type == "Promoter") %>%
  mutate(Peak_Type = "Enhancer") %>%
  select(peak, cicero_Gene, Peak_Type) %>%
  rename(Peak = peak) %>%
  rename(Gene = cicero_Gene) %>%
  mutate(Is.CAP = T) %>%
  mutate(Is.DAP = T)
  
# get unique co-accessible peaks--------------------------------------------------------------
coaccess_uniq_df <- coaccess_filtered_df %>%
  filter(Type == "Promoter" | cicero_Type == "Promoter") %>%
  select(peak_coacess) %>%
  unique()
peak2gene_df <- rbind(coaccess_enhancer_df, coaccess_promoter_df, dap_enhancer_df)
peak2gene_df <- unique(peak2gene_df)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Coaccessible_Peaks.", run_id, ".tsv")
write.table(x = coaccess_uniq_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Peak2Gene.", run_id, ".tsv")
write.table(x = peak2gene_df, file = file2write, quote = F, sep = "\t", row.names = F)
