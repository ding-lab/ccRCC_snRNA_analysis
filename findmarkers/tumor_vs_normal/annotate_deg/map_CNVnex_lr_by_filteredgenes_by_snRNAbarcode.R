# Yige Wu @ WashU 2021 Jun
## annotate sample copy number profile
## CNV: https://wustl.box.com/s/vlde6w791k81q1aibvgs0a487zxleij6
# Purity and ploidy: https://wustl.box.com/s/jel5krgvnvlq5z32vdg4itdcq6znuqzn
# From UMich

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
## input CNA matrix
cna_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Combined/Absolute_cnv/c3-ccrcc-combined-cnvex-lr_v1.0.csv")
## input snRNA sample set
metadata_df <- fread("./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv", data.table = F)
## input the barcode info
barcodes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv")
## input prefiltered genes to test
genes_filtered_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/Tumorcells_vs_PTcells/20210823.v1/Each31RCC_Tumorcells_vs_CombinedPTcells.20210823.v1.tsv")

# preprocess ----------------------------
## get barcodes to process
barcodes_df <- barcodes_df %>%
  mutate(id_bc = paste0(orig.ident, "_", individual_barcode))
barcodes_df$Case <- mapvalues(x = barcodes_df$orig.ident, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Case))
barcodes_df$Sample_type <- mapvalues(x = barcodes_df$orig.ident, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Sample_Type))
barcodes_df <- barcodes_df %>%
  mutate(keep = ((Sample_type == "Tumor" & Cell_group_w_epithelialcelltypes == "Tumor cells") | (Sample_type == "Normal" & Cell_group_w_epithelialcelltypes == "Proximal tubule")))
table(barcodes_df$Cell_type.detailed[barcodes_df$keep])
# barcodes_df %>%
#   filter(Sample_type == "Tumor" & Cell_type == "PT")
barcodes_process <- barcodes_df$id_bc[barcodes_df$keep]
cases_process <- barcodes_df$Case[barcodes_df$keep]
celltypes_process <- barcodes_df$Cell_group_w_epithelialcelltypes[barcodes_df$keep]
table(celltypes_process)
## get genes to process
genes_filtered_df <- genes_filtered_df %>%
  filter(easyid_tumor != "C3L-00359-T1")
genes_process <- unique(genes_filtered_df$genesymbol_deg)
genes_process <- genes_process[genes_process %in% cna_df$gene_name]
length(genes_process)

# preprocess mean CNV values per gene--------------------------------------------------------------
## preprocess the CNV data frame
colnames_old <- colnames(cna_df)
colnames_new <- str_split_fixed(string = colnames_old, pattern = "\\.", n = 4)[,1]
colnames(cna_df) <- colnames_new
## filter the CNVs
cna_filtered_df <- cna_df %>%
  filter(gene_name %in% genes_process)
cna_filtered_df <- cna_filtered_df[!duplicated(cna_filtered_df$gene_name),]
cna_bybc_df <- cna_filtered_df[, cases_process]
cna_bybc_df[, celltypes_process == "PT"] <- 0
rownames(cna_bybc_df) <- cna_filtered_df$gene_name
colnames(cna_bybc_df) <- barcodes_process
cna_bybc_df[1:5, 1:5]
cna_t_mat <- t(as.matrix(cna_bybc_df))
cna_t_mat[1:5, 1:5]

# write table -------------------------------------------------------------
file2write <- paste0(dir_out, "Barcode2Gene.CNV.", run_id, ".RDS")
saveRDS(object = cna_t_mat, file = file2write, compress = T)
