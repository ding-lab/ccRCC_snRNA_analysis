# Yige Wu @WashU March 2020
## for calculating the fraction of tumor cells with cnv in different frequently altered genes

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode to chr-region cnv state
barcode2cnv_bychr_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/annotate_barcode_cnv/20200228.v1/Expected_CNV_State_By_Chr_Region_By_Barcode.20200228.v1.tsv", data.table = F)
## input barcode to cell type info
barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/map_celltype_to_barcode/20200309.v1/30_aliquot_integration.barcode2celltype.20200309.v1.tsv", data.table = F)
## input id meta data
id_metadata_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## chromosome regions from which the CNVs are annotated here
chr_regions2process <- unique(ccrcc_cna_genes_df$chr_region)
chr_regions2process <- as.vector(chr_regions2process)

# filter to only tumor cells and only tumor sample ----------------------------------------------
## add integrated barcode to the cnv table
barcode2cnv_bychr_df <- merge(barcode2cnv_bychr_df %>%
                                rename(individual_barcode = barcode), 
                              barcode2celltype_df %>%
                                select(orig.ident, integrated_barcode, individual_barcode) %>%
                                rename(aliquot = orig.ident),
                              by = c("aliquot", "individual_barcode"),
                              all.x = T)
## get the tumor aliquot ids
tumor_aliquot_ids <- id_metadata_df$Aliquot.snRNA[id_metadata_df$Sample_Type == "Tumor"]
## get malignant nephron epithelium cell barcodes
malignant_barcodes <- barcode2celltype_df$integrated_barcode[barcode2celltype_df$Most_Enriched_Cell_Group == "Nephron_Epithelium" & barcode2celltype_df$Is_Normal_Nephron_Epithelium == F]
## filter by malignant cell barcodes and tumor samples only
tumorcell2cnv_bychr_df <- barcode2cnv_bychr_df %>%
  filter(integrated_barcode %in% malignant_barcodes) %>%
  filter(aliquot %in% tumor_aliquot_ids)
# calculate the fraction of tumor cells with expected cnv -----------------
## group by sample and aggregate the number of barcodes with expected cnv
tumorcell2cnv_bychr_mat <- tumorcell2cnv_bychr_df[, chr_regions2process]
tumorcell2isexpectedcnv_bychr_mat <- (!is.na(tumorcell2cnv_bychr_mat) & tumorcell2cnv_bychr_mat == "Expected") + 0
count_expectedcnv_bychr_byaliquot_df <- aggregate(x = tumorcell2isexpectedcnv_bychr_mat, list(tumorcell2cnv_bychr_df$aliquot), sum)
count_expectedcnv_bychr_byaliquot_df <- count_expectedcnv_bychr_byaliquot_df %>%
  rename(aliquot = Group.1)
## group by sample and aggregate the number of total barcodes
tumorcell2iscnv_bychr_mat <- (!is.na(tumorcell2cnv_bychr_mat)) + 0
count_cnv_bychr_byaliquot_df <- aggregate(x = tumorcell2iscnv_bychr_mat, list(tumorcell2cnv_bychr_df$aliquot), sum)
count_cnv_bychr_byaliquot_df <- count_cnv_bychr_byaliquot_df %>%
  rename(aliquot = Group.1)
# write outputs -------------------------------------------
## write the number cells with expected cnv
file2write <- paste0(dir_out, "number_of_tumorcells.expectedCNA.by_chr_region.", run_id, ".tsv")
write.table(x = count_expectedcnv_bychr_byaliquot_df, file = file2write, quote = F, sep = "\t", row.names = F)
## write the number cells with cnv calls
file2write <- paste0(dir_out, "number_of_tumorcells.withCNAcalls.by_chr_region.", run_id, ".tsv")
write.table(x = count_cnv_bychr_byaliquot_df, file = file2write, quote = F, sep = "\t", row.names = F)
## write the number cells with cnv calls
frac_expectedcnv_bychr_byaliquot_df <- count_expectedcnv_bychr_byaliquot_df[, chr_regions2process]/count_cnv_bychr_byaliquot_df[,chr_regions2process]
frac_expectedcnv_bychr_byaliquot_df <- cbind(data.frame(aliquot = count_expectedcnv_bychr_byaliquot_df$aliquot),
                                             frac_expectedcnv_bychr_byaliquot_df)
file2write <- paste0(dir_out, "fraction_of_tumorcells.expectedCNA.by_chr_region.", run_id, ".tsv")
write.table(x = frac_expectedcnv_bychr_byaliquot_df, file = file2write, quote = F, sep = "\t", row.names = F)
