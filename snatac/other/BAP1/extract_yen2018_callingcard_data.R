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
## input the calling card data
bap1_pbase_df <- readxl::read_excel(path = "./Literature_Review/Genes/BAP1/2018_BMC_BAP1_Uveal_Melanoma/12920_2018_424_MOESM1_ESM.xlsx")
pbase_bap1_df <- readxl::read_excel(path = "./Literature_Review/Genes/BAP1/2018_BMC_BAP1_Uveal_Melanoma/12920_2018_424_MOESM2_ESM.xlsx")
## input the our bap1 hits
bap1_tophits_df <-  readxl::read_excel(path = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Peak_Annotation/BAP1/BAP1_DACR_genes_overlap_m1_probe_genes.xlsx")

# filter ------------------------------------------------------------------
peaks_df <- rbind(bap1_pbase_df %>%
                    mutate(construct = "BAP1-PBase"),
                  pbase_bap1_df %>%
                    mutate(construct = "PBase-BAP1"))
chromatin_filtered_df <- peaks_df %>%
  filter(`Nearest Feature Name` %in% bap1_tophits_df$gene_symbol)
bs_long_df <- chromatin_filtered_df %>%
  group_by(`Nearest Feature Name`, construct) %>%
  summarise(num_bindingsites = n()) %>%
  rename(gene_symbol.peak = `Nearest Feature Name`)
bs_wide_df <- dcast(data = bs_long_df, formula = gene_symbol.peak ~ construct)
rownames(bs_wide_df) <- bs_wide_df$gene_symbol.peak
bs2merge_df <- bs_wide_df[bap1_tophits_df$gene_symbol,]

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "2018_BMC_BAP1_Uveal_Melanoma.", "bindingsites", ".tsv")
write.table(file = file2write, x = bs2merge_df, quote = F, sep = "\t", row.names = F)
