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
deg_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/Tumor_Subcluster_Group/20210421.v1/Immune.TumorManualCluster.DEGs.Wilcox.Minpct0.1.Logfc0.min.diff.pct0.1.tsv")
## input pathway to genes
genes_immune_df <- fread(data.table = F, input = "./Resources/Knowledge/Databases/ImmPort/GeneList.txt")

# overlap -----------------------------------------------------------------
deg_sig_df <- deg_df %>%
  filter(p_val_adj < 0.05) %>%
  filter(avg_logFC > 0) %>%
  filter(genesymbol_deg %in% genes_immune_df$Symbol)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ImmuneModuleEnrichedDEGs.Filtered.tsv")
write.table(x = deg_sig_df, file = file2write, quote = F, sep = "\t", row.names = F)

