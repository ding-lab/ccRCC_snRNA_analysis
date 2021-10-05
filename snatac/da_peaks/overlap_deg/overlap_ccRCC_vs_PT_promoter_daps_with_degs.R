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
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_ccrcc_vs_pt_promoter_daps/20211004.v1/ccRCC_vs_PTDAPs.Annotated.Promoter.20211004.v1.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_individual_and_CNVcorrected_DEGs/20210824.v1/Consistent.Tumor_vs_PT_DEGs.CNVcorrected.20210824.v1.tsv")

# annotate and filter peaks ------------------------------------------------------------
peaks2degs_df <- merge(x = peaks_anno_df, 
                       y = degs_df %>%
                         rename(avg_log2FC = avg_log2FC.allTumorcellsvsPT) %>%
                         dplyr::select(genesymbol_deg, Num_sig_up, Num_sig_down, avg_log2FC),
                       by.x = c("Gene"), by.y = c("genesymbol_deg"), suffix = c(".snATAC", ".snRNA"))
peaks2degs_filtered_df <- peaks2degs_df %>%
  dplyr::filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA))

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_vs_PT_DEG_associated_DAP2DEG.", run_id, ".tsv")
write.table(x = peaks2degs_df, file = file2write, quote = F, sep = "\t", row.names = F)

peaks2degs_df %>%
  filter(peak2gene_type == "Promoter") %>%
  filter(DAP_direction == "Up") %>%
  filter(avg_log2FC.snRNA > 0) %>%
  select(peak) %>%
  unique() %>%
  nrow()

peaks2degs_df %>%
  filter(peak2gene_type == "Promoter") %>%
  filter(DAP_direction == "Down") %>%
  filter(avg_log2FC.snRNA < 0) %>%
  select(peak) %>%
  unique() %>%
  nrow()

test_df <- peaks2degs_df %>%
  filter(peak2gene_type == "Promoter") %>%
  filter(DAP_direction == "Up") %>%
  filter(avg_log2FC.snRNA > 0)

peaks2degs_df %>%
  filter(peak2gene_type == "Enhancer") %>%
  filter(DAP_direction == "Up") %>%
  filter(avg_log2FC.snRNA > 0) %>%
  select(peak) %>%
  unique() %>%
  nrow()
