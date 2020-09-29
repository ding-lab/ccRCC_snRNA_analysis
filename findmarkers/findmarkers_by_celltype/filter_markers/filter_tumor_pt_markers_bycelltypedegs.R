# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the marker genes
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/findallmarker_wilcox_cellgroup_on_katmai/20200714.v2/findallmarkers_wilcox_bycellgroup.pos..logfcthreshold0.1.minpct0.1.mindiffpct0.1.tsv")
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200911.v1.tsv")

# get top 50 genes by avg_logFC ------------------------------------------------
deg_filtered_df <- deg_df %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(cluster == "Tumor cells") %>%
  filter(gene %in% gene2celltype_df$Gene[gene2celltype_df$Cell_Type1 %in% c("Tumor cells", "Proximal tubule", "Fibroblasts")])
deg_filtered_df$Cell_Type1 <- mapvalues(x = deg_filtered_df$gene, from = gene2celltype_df$Gene, to = as.vector(gene2celltype_df$Cell_Type1))
deg_filtered_df <- deg_filtered_df %>%
  group_by(Cell_Type1) %>%
  top_n(n = 3, wt = avg_logFC)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Essential_Tumor_Cell_Markers", ".tsv")
write.table(x = deg_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

