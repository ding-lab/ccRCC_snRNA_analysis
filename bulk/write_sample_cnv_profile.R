# Yige Wu @ WashU 2020 Feb
## annotate sample copy number profile (3p, 5q, 14q)

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
## input CNA matrix
cna_tab <- loadCNAstatus()
## input snRNA sample set
# snRNA_meta_data <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_cluster/write_path_to_seurat_objects_on_box/20200219.v1/Seurat_Object_Paths.20200219.v1.tsv", data.table = F)
## specify the chromosome regions to process
chr_regions <- c("3p", "5q", "14q")
# select genes to evaluate copy number profile ----------------------------
## get the aliquot ids to process
aliquots <- unique(srat_paths$Case)
## specify the genes to be used 
genes2plot <- ccrcc_cna_genes_df$gene_symbol[ccrcc_cna_genes_df$chr_region %in% c("3p", "5q", "14q")]
## filter the CNVs
cna_tab <- cna_tab[cna_tab$gene %in% genes2plot, c("gene", aliquots)]
rownames(cna_tab) <- cna_tab$gene

# for each chromosome region, get the cnv state for each aliquot ----------
chr_cnv_state_list <- lapply(chr_regions, function(chr_region, cnv_df, gene2chr_df) {
  genes_tmp <- as.vector(gene2chr_df$gene_symbol[gene2chr_df$chr_region == chr_region])
  loss_count <- colSums(cnv_df[genes_tmp,-1] == "deletion")
  gain_count <- colSums(cnv_df[genes_tmp,-1] == "amplification")
  chr_cnv_state <- rep("loss", length(loss_count))
  chr_cnv_state[gain_count == 0 & loss_count == 0] <- ""
  chr_cnv_state[gain_count > 0 & loss_count == 0] <- "gain"
  chr_cnv_state[gain_count > 0 & loss_count > 0] <- "mixed"
  return(chr_cnv_state)
}, cnv_df = cna_tab, gene2chr_df = ccrcc_cna_genes_df)
## make the list into matrix
chr_cnv_state_mat <- matrix(data = unlist(chr_cnv_state_list), ncol = length(chr_regions), 
                                dimnames = list(row_names = colnames(cna_tab[,-1]),
                                                col_names = chr_regions))
chr_cnv_state_mat
## make matrix into data frame
chr_cnv_state_df <- data.frame(Case = rownames(chr_cnv_state_mat))
chr_cnv_state_df <- cbind(chr_cnv_state_df, as.data.frame(chr_cnv_state_mat))
chr_cnv_state_df$chr_cnvs_summary <- paste0(ifelse(chr_cnv_state_df$`3p` != "", paste0("|", "3p_", chr_cnv_state_df$`3p`, "|"), ""),
                                            ifelse(chr_cnv_state_df$`5q` != "", paste0("|","5q_", chr_cnv_state_df$`5q`, "|"), ""),
                                            ifelse(chr_cnv_state_df$`14q` != "", paste0("|", "14q_", chr_cnv_state_df$`14q`, "|"), ""))


# write table -------------------------------------------------------------
write.table(x = chr_cnv_state_df, file = paste0(dir_out, "Bulk_WGS_Chr_CNV_Profile", ".", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

