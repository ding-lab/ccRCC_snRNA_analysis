# Yige Wu @WashU Oct 2020

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
## input the DEG-TF matrix
deg2tf_wide_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/DEGs_with_TFs_inDARs.tsv")
## input the deg table
sn_rna_de_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_tumor_vs_pt_on_katmai/20200929.v1/findallmarkers_wilcox_tumorcells_vs_pt.20200929.v1.tsv")
## input the differentially expression of bulk 
bulk_rna_de_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/Tumor_normal_bulk_RNA.txt")
bulk_pro_de_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/Tumor_normal_bulk_Protein.txt")

# merge with snRNA DEGs ---------------------------------------------------
deg2tf_de_df <- deg2tf_wide_df %>%
  rename(genesymbol_deg = V1)
## merge with snRNA
colnames(sn_rna_de_df) <- paste0(colnames(sn_rna_de_df), ".snrna")
deg2tf_de_df <- merge(x = deg2tf_de_df, y = sn_rna_de_df, 
                      by.x = c("genesymbol_deg"), by.y = c("row_name.snrna"), all.x = T)

## merge with bulk mRNA
colnames(bulk_rna_de_df) <- paste0(colnames(bulk_rna_de_df), ".bulk.rna")
deg2tf_de_df <- merge(x = deg2tf_de_df, y = bulk_rna_de_df, 
                      by.x = c("genesymbol_deg"), by.y = c("Gene.bulk.rna"), all.x = T)
## merge with bulk protein
colnames(bulk_pro_de_df) <- paste0(colnames(bulk_pro_de_df), ".bulk.pro")
deg2tf_de_df <- merge(x = deg2tf_de_df, y = bulk_pro_de_df, 
                      by.x = c("genesymbol_deg"), by.y = c("Gene.bulk.pro"), all.x = T)
## sort
deg2tf_de_df <- deg2tf_de_df %>%
  arrange(FDR.bulk.pro)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs_with_TFs_inDARs", ".sn_and_Bulk_DE_Annotated", ".tsv")
write.table(x = deg2tf_de_df, file = file2write, quote = F, sep = "\t", row.names = F)
