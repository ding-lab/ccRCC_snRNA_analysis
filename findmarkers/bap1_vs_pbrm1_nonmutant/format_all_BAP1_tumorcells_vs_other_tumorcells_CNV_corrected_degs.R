# Yige Wu @WashU Jun 2021

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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/findallmarker_LR_all_BAP1_tumorcells_vs_PBRM1_NonMutant_cells_on_katmai/20210609.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")
deg_cnvcor_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/findmarker_LR_wCNV_all_BAP1_tumorcells_vs_PBRM1_NonMutant_cells_on_katmai/20210609.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")

# merge -------------------------------------------------------------------
deg_new_df <- merge(x = deg_cnvcor_df %>%
                      rename(genesymbol_deg = row.names), 
                    y = deg_df %>%
                      select(genesymbol_deg, avg_log2FC, pct.1, pct.2),
                    by = c("genesymbol_deg"), all.x = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.CNVcorrected.tsv")
write.table(x = deg_new_df, file = file2write, quote = F, sep = "\t", row.names = F)

deg_new_df %>%
  filter(FDR < 0.05) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down")) %>%
  select(direction) %>%
  table()
