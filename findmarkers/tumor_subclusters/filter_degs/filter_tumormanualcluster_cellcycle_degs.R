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
deg_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/Tumor_Subcluster_Group/20210421.v1/Cell_cycle.TumorManualCluster.DEGs.Wilcox.Minpct0.1.Logfc0.min.diff.pct0.1.tsv")
## input pathway to genes
pathway2genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/map_genes2pathway_w_avgexp_sct_data/20210414.v1/TumorCluster.DEG2Pathway.tsv")

# overlap -----------------------------------------------------------------
deg_sig_df <- deg_df %>%
  filter(p_val_adj < 0.05) %>%
  filter(avg_logFC > 0)

deg_int_df <- merge(x = deg_sig_df, y = pathway2genes_df, by.x = c("genesymbol_deg"), by.y = c("GeneSymbol"))
deg_int_df <- deg_int_df %>%
  arrange(GeneSet_Name, desc(avg_logFC))
## filter to specific pathways
deg_int_filtered_df <- deg_int_df %>%
  filter(GeneSet_Name %in% c("HALLMARK_MITOTIC_SPINDLE", "HALLMARK_G2M_CHECKPOINT")) %>%
  arrange(desc(GeneSet_Name))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellCycleModuleEnrichedDEGs.Filtered.tsv")
write.table(x = deg_int_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

