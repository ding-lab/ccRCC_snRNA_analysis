# Yige Wu @WashU March 2020
## for plotting the marker genes for integrated object, showing cell of origin

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
## input the barcode to CNV state info
barcode2cnv_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/annotate_barcode_cnv/20200228.v1/Expected_CNV_State_By_Chr_Region_By_Barcode.20200228.v1.tsv", data.table = F)
## input barcode to cell type info
barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/map_celltype_to_barcode/20200224.v1/30_aliquot_integration.barcode2celltype.20200224.v1.tsv", data.table = F)

# plot by each aliquot ----------------------------------------------------
aliquot_tmp <- "CPT0001220012"

### get malignant nephron epithelium cell barcodes
malignant_barcodes <- barcode2celltype_df$individual_barcode[barcode2celltype_df$orig.ident == aliquot_tmp & barcode2celltype_df$Most_Enriched_Cell_Group == "Nephron_Epithelium" & barcode2celltype_df$Is_Normal_Nephron_Epithelium == F]

## filter the cells by aliquot
barcode2cnv_aliquot_df <- barcode2cnv_df %>%
  filter(aliquot == aliquot_tmp) %>%
  filter(barcode %in% malignant_barcodes)

## create the united cnv id to make the counting easier
barcode2cnv_aliquot_df$united_cnv_id <- sapply(1:nrow(barcode2cnv_aliquot_df), function(i, input_data_frame) {
  text <- paste0(input_data_frame[i,], collapse = "|")
  return(text)
}, input_data_frame = barcode2cnv_aliquot_df[,as.vector(unique(ccrcc_cna_genes_df$chr_region))])

## count the occurence of the united cnv id
cnv_cooccurrence_barcode_count <- barcode2cnv_aliquot_df %>%
  select(united_cnv_id) %>%
  table() %>%
  as.data.frame() %>%
  arrange(desc(Freq))

## recover the CNV pattern
### rename the column names
colnames(cnv_cooccurrence_barcode_count) <- c("united_cnv_id", "Freq")
### get the corespinding CNV pattern for each united cnv id
united_cnv_id_df <- str_split_fixed(string = cnv_cooccurrence_barcode_count$united_cnv_id, pattern = "\\|", n = length(as.vector(unique(ccrcc_cna_genes_df$chr_region))))
united_cnv_id_df <- as.data.frame(united_cnv_id_df)
colnames(united_cnv_id_df) <- as.vector(unique(ccrcc_cna_genes_df$chr_region))
### bind with the frequency count
cnv_cooccurrence_barcode_count <- cbind(cnv_cooccurrence_barcode_count, united_cnv_id_df)


barcode2cnv_aliquot_df.m <- melt(barcode2cnv_aliquot_df %>%
                                   select(-united_cnv_id), id.vars = c("aliquot", "barcode"))
cnv_count <- barcode2cnv_aliquot_df.m %>%
  filter(value == "Expected") %>%
  select(variable, value) %>%
  table() %>%
  as.data.frame()
