# Yige Wu @ WashU 2020 Feb
## annotate sample copy number profile (3p, 5q, 14q)
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
metadata_df <- fread("./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv", data.table = F)
## input ccRCC CNA genes
ccrcc_cna_genes_df <- readxl::read_excel(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## specify the chromosome regions to process
chr_regions <- c("3p", "5q", "14q", "7p")
# select genes to evaluate copy number profile ----------------------------
## get the aliquot ids to process
aliquots <- unique(metadata_df$Case[metadata_df$snRNA_available])

# preprocess --------------------------------------------------------------
ccrcc_cna_genes_df <- ccrcc_cna_genes_df %>%
  mutate(Chromosome_Arm = paste0(Chromosome, Arm))
## specify the genes to be used 
genes2plot <- ccrcc_cna_genes_df$Gene_Symbol[ccrcc_cna_genes_df$Chromosome_Arm %in% chr_regions]
## filter the CNVs
cna_filtered_df <- cna_df[cna_df$gene_name %in% genes2plot,]
colnames_old <- colnames(cna_filtered_df)
colnames_new <- str_split_fixed(string = colnames_old, pattern = "\\.", n = 4)[,1]
colnames(cna_filtered_df) <- colnames_new
rownames(cna_filtered_df) <- cna_filtered_df$gene_name

# for each chromosome region, get the cnv state for each aliquot ----------
chr_cnv_state_list <- lapply(chr_regions, function(chr_region, cnv_df, gene2chr_df) {
  genes_tmp <- as.vector(gene2chr_df$Gene_Symbol[gene2chr_df$Chromosome_Arm == chr_region])
  cnv_value <- colMeans(cnv_df[genes_tmp, 3:ncol(cnv_df)])
  return(cnv_value)
}, cnv_df = cna_filtered_df, gene2chr_df = ccrcc_cna_genes_df)
## make the list into matrix
chr_cnv_state_mat <- matrix(data = unlist(chr_cnv_state_list), ncol = length(chr_regions), 
                                dimnames = list(row_names = colnames(cna_filtered_df[,3:ncol(cna_filtered_df)]),
                                                col_names = chr_regions))
View(chr_cnv_state_mat)
colnames(chr_cnv_state_mat) <- paste0("CN.", colnames(chr_cnv_state_mat))
## make matrix into data frame
chr_cnv_state_df <- data.frame(Case = rownames(chr_cnv_state_mat))
chr_cnv_state_df <- cbind(chr_cnv_state_df, as.data.frame(chr_cnv_state_mat))

# write table -------------------------------------------------------------
write.table(x = chr_cnv_state_df, file = paste0(dir_out, "Bulk_WGS_Chr_CNV_Values.", "UsingKnownCNVGenes.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

