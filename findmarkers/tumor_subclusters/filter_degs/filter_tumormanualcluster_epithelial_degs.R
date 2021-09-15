# Yige Wu @WashU Apr 2021

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
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_tumor_w_enriched_genesetmodule_vs_others_doparallel_katmai/20210624.v25Cores/EMT.TumorManualCluster.DEGs.Wilcox.Minpct0.1.Logfc0.min.diff.pct0.1.tsv")
## input known EMT markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Human.Gene2CellType.20210819.tsv")

# overlap -----------------------------------------------------------------
deg_sig_df <- deg_df %>%
  filter(p_val_adj < 0.05)

deg_epithelia_df <- deg_sig_df %>%
  filter(avg_log2FC < 0 ) %>%
  filter(genesymbol_deg %in% emt_genes_df$Gene[emt_genes_df$Gene_Group2 %in% c("Epithelial", "Proximal tubule")]) %>%
  arrange(avg_log2FC, p_val_adj)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellCycleModuleEnrichedDEGs.Filtered.tsv")
write.table(x = deg_int_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

