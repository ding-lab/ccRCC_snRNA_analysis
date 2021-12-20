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
coaccess_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Coaccessibility/28_ccRCC_snATAC_CICERO.0.25_cutoff.tsv")
## input peak2gene annotation
peak2gene_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Peak_Annotation/28_snATACmerged_allPeaks.Annotated.20210712.tsv")

# filter and extract ------------------------------------------------------
nrow(coaccess_df)
coaccess_filtered_df <- coaccess_df %>%
  filter(coaccess >= 0.25)
nrow(coaccess_filtered_df)
## prepare gene annotation
peak2gene_df <- peak2gene_df %>%
  mutate(peak2gene_type = str_split_fixed(string = Type, pattern = " \\(", n = 2)[,1])
table(peak2gene_df$peak2gene_type)
### filter to only peak 1
peak2gene_df1 <- peak2gene_df %>%
  filter(peak %in% coaccess_filtered_df$Peak1) %>%
  mutate(genesymbol = Gene) %>%
  select(peak, peak2gene_type, genesymbol)
peak2gene_df2 <- peak2gene_df %>%
  filter(peak %in% coaccess_filtered_df$Peak2) %>%
  mutate(genesymbol = Gene) %>%
  select(peak, peak2gene_type, genesymbol)
## merge
coaccess_peaks2genes_df <- merge(x = coaccess_filtered_df, y = peak2gene_df1, by.x = c("Peak1"), by.y = c("peak"), all.x = T)
coaccess_peaks2genes_df <- merge(x = coaccess_peaks2genes_df, y = peak2gene_df2, by.x = c("Peak2"), by.y = c("peak"), all.x = T, suffixes = c(".1", ".2"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Coaccessible_Peaks.Annotated.", run_id, ".tsv")
write.table(x = coaccess_peaks2genes_df, file = file2write, quote = F, sep = "\t", row.names = F)
