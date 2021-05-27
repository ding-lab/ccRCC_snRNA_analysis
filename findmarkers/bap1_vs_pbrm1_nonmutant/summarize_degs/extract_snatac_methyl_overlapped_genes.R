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
exp_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/summarize_BAP1_vs_PBRM1_NonMutant_DEGs/20210430.v3/BAP1_DEGs.Sig.20210430.v3.tsv")
## input the our bap1 hits
bap1_tophits_df <-  readxl::read_excel(path = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Peak_Annotation/BAP1/BAP1_DACR_genes_overlap_m1_probe_genes.xlsx")

# filter ------------------------------------------------------------------
exp_filtered_df <- exp_df %>%
  filter(genesymbol_deg %in% bap1_tophits_df$gene_symbol) %>%
  unique()
exp_filtered_df <- as.data.frame(exp_filtered_df)
rownames(exp_filtered_df) <- exp_filtered_df$genesymbol_deg
exp2merge_df <- exp_filtered_df[bap1_tophits_df$gene_symbol,]

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "snRNA.", ".tsv")
write.table(file = file2write, x = exp2merge_df, quote = F, sep = "\t", row.names = F)
