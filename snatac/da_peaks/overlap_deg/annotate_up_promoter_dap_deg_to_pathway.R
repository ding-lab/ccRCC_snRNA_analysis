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
##input overlap
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/overlap_deg/overlap_up_dap_w_up_snRNA_bulk_degs/20210510.v1/Up_Promoter_DAP.Overlap.Up_snRNA_bulk_DEGs.20210510.v1.tsv")
## input pathway annotation
gene2pathway_df1 <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/map_genes2pathway__tumorcells_up_degs_msigdb_H_CP_sig_pathways/20210430.v1/DEG2Pathway.Collapsed.tsv")
gene2pathway_df2 <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/map_genes2pathway_tumorcells_up_degs_msigdb_H_CP/20210430.v1/DEG2Pathway.Collapsed.tsv")

# merge -------------------------------------------------------------------
deg_anno_df <- merge(x = deg_df, y = gene2pathway_df1, by.x = c("genesymbol_deg"), by.y = c("GeneSymbol"), all.x = T)
deg_anno_df <- merge(x = deg_anno_df, y = gene2pathway_df2, by.x = c("genesymbol_deg"), by.y = c("GeneSymbol"), all.x = T, suffixes = c("1", "2"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Up_Promoter_DAP.Overlap.Up_snRNA_bulk_DEGs.20210510.v1.tsv")
write.table(x = deg_anno_df, file = file2write, quote = F, sep = "\t", row.names = F)
