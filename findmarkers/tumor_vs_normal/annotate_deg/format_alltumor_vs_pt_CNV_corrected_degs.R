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
deg_cnvcor_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_LR_wCNV_all_ccRCC_vs_pt_on_katmai/20210823.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.1.AssayRNA.tsv")
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_LR_all_ccRCC_vs_pt_on_katmai/20210823.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.1.AssayRNA.tsv")
cnv_per_feature_df=readRDS('./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/map_CNVnex_lr_by_filteredgenes_by_snRNAbarcode_katmai/20210823.v1/Barcode2Gene.CNV.20210823.v1.RDS')

# process -----------------------------------------------------------------
colnames(deg_cnvcor_df) <- paste0(colnames(deg_cnvcor_df), ".CNVcorrected")
deg_cnvcor_df$genesymbol_deg <- colnames(cnv_per_feature_df)
rm(cnv_per_feature_df)
colnames(deg_df) <- paste0(colnames(deg_df), ".allTumorcellsvsPT")
deg_merged_df <- merge(x = deg_cnvcor_df, 
                       y = deg_df %>%
                         select(genesymbol_deg.allTumorcellsvsPT, avg_log2FC.allTumorcellsvsPT, pct.1.allTumorcellsvsPT, pct.2.allTumorcellsvsPT), 
                       by.x = c("genesymbol_deg"), by.y = c("genesymbol_deg.allTumorcellsvsPT"), all.x = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "LR.logfc.threshold0.min.pct0.1.min.diff.pct0.1.AssayRNA.tsv")
write.table(x = deg_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
