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
dap2gene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_peaks_with_enhancer_prediction/20210514.v1/DAPeak2Gene.Count5.20210514.v1.tsv")
## input coaccessiblity peaks
cap2gene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/get_unique_peaks_from_coaccessibility_analysis/20210517.v1/Peak2Gene.20210517.v1.tsv")
## input known enhancers
enhancers_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/Knowledge/EnhancerAtlast/species_enh_csv/matrix_hs.csv")

# merge -------------------------------------------------------------------
peak2gene_df <- rbind(dap2gene_df %>%
                        mutate(Is.CAP = (Peak %in% cap2gene_df$Peak[cap2gene_df$Is.CAP])) %>%
                        mutate(Is.DAP = T), 
                      cap2gene_df)
peak2gene_df <- unique(peak2gene_df)
table(peak2gene_df$Is.DAP)
## some 
peak2gene_df$Is.DAP[peak2gene_df$Peak %in% dap2gene_df$Peak] <- T
peak2gene_df <- unique(peak2gene_df)

table(peak2gene_df$Is.DAP)
peak2gene_df %>%
  select(Is.DAP, Is.CAP) %>%
  table()
peak_uniq_df <- peak2gene_df %>%
  select(Peak) %>%
  unique()
nrow(peak_uniq_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Peaks.", run_id, ".tsv")
write.table(x = peak_uniq_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Peak2Gene.", run_id, ".tsv")
write.table(x = peak2gene_df, file = file2write, quote = F, sep = "\t", row.names = F)
