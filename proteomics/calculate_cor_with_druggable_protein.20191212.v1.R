# Yige Wu @ WashU 2019 Dec
## do global correlation for specific druggable gene expression/protein/phosphoprotein

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# specify the least number of datapoints for correlation ------------------
least_no_datapoints <- 20

# specify gene and expression level to do correlation with all other genes --------
geneA <- "VEGFA"
geneA_exp_type <- "PRO"


# input bulk meta data ----------------------------------------------------
bulk_meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")

# get aliquot ids for tumor samples ---------------------------------------
tumor_bulk_aliquot_ids <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Tumor"]
tumor_bulk_aliquot_ids

# input bulk protein data ------------------------------------------------------
protein_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/proteome/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)

# input bulk phosphosite data ------------------------------------------------------
phospho_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/phosphoproteome/6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv", data.table = F)

# input RNA data ----------------------------------------------------------
rna_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/mRNA/RNA_rpkm_tumor_normal.tsv", data.table = F)
tumor_bulk_rna_ids <- bulk_meta_tab$RNA.ID[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Tumor"]
rna_tab <- rna_tab[, c("geneID", tumor_bulk_rna_ids)]
colnames(rna_tab) <- c("geneID", mapvalues(x = tumor_bulk_rna_ids, from = bulk_meta_tab$RNA.ID, to = as.vector(bulk_meta_tab$Specimen.Label)))

# get expression level of the gene to test --------------------------------
if (geneA_exp_type == "PRO") {
  geneA_exp_values <- protein_tab[protein_tab$Index == geneA, tumor_bulk_aliquot_ids]
}

geneA_exp_nonna_ids <- colnames(geneA_exp_values)[!is.na(geneA_exp_values)]

# test correlation with mRNA levels of other genes ---------------------
row_indexs <- which(rowSums(rna_tab[,geneA_exp_nonna_ids] > 0) >= least_no_datapoints)

sample_sizes <- vector(mode = "numeric", length = length(row_indexs))
p_values <- vector(mode = "numeric", length = length(row_indexs))
rhos <- vector(mode = "numeric", length = length(row_indexs))

for (i in 1:length(row_indexs)) {
  geneB <- rna_tab[row_indexs[i], "geneID"]
  geneB_exp_values <- rna_tab[row_indexs[i], tumor_bulk_aliquot_ids]
  
  nonna_ids <- tumor_bulk_aliquot_ids[!is.na(geneA_exp_values) & geneB_exp_values > 0]
  
  cor_stats_tmp <- cor.test(x = unlist(geneA_exp_values[,nonna_ids]), y = unlist(geneB_exp_values[,nonna_ids]), method = "spearman")
  p_values[i] <- cor_stats_tmp$p.value
  rhos[i] <- cor_stats_tmp$estimate
  sample_sizes[i] <- length(nonna_ids)
  print(i)
}

cor_stats_rna_tab <- data.frame(geneB = rna_tab[row_indexs, "geneID"],
                                geneB_exp_type = "RNA", 
                                geneB_exp_site = "",
                                genepair_cor_p_value = p_values,
                                genepair_cor_rho = rhos,
                                genepair_cor_sample_size = sample_sizes,
                                geneA = geneA, 
                                geneA_exp_type = geneA_exp_type)
cor_stats_rna_tab$genepair_cor_fdr <- p.adjust(p = cor_stats_rna_tab$genepair_cor_p_value, method = "fdr")


# test correlation with protein levels of other genes ---------------------
row_indexs <- which(rowSums(!is.na(protein_tab[,geneA_exp_nonna_ids])) >= least_no_datapoints)

sample_sizes <- vector(mode = "numeric", length = length(row_indexs))
p_values <- vector(mode = "numeric", length = length(row_indexs))
rhos <- vector(mode = "numeric", length = length(row_indexs))

for (i in 1:length(row_indexs)) {
  geneB <- protein_tab[row_indexs[i], "Index"]
  geneB_exp_values <- protein_tab[row_indexs[i], tumor_bulk_aliquot_ids]
  
  nonna_ids <- tumor_bulk_aliquot_ids[!is.na(geneA_exp_values) & !is.na(geneB_exp_values)]
  
  cor_stats_tmp <- cor.test(x = unlist(geneA_exp_values[,nonna_ids]), y = unlist(geneB_exp_values[,nonna_ids]), method = "spearman")
  p_values[i] <- cor_stats_tmp$p.value
  rhos[i] <- cor_stats_tmp$estimate
  sample_sizes[i] <- length(nonna_ids)
  print(i)
}

cor_stats_pro_tab <- data.frame(geneB = protein_tab[row_indexs, "Index"],
                                geneB_exp_type = "PRO", 
                                geneB_exp_site = "",
                                genepair_cor_p_value = p_values,
                                genepair_cor_rho = rhos,
                                genepair_cor_sample_size = sample_sizes,
                                geneA = geneA, 
                                geneA_exp_type = geneA_exp_type)
cor_stats_pro_tab$genepair_cor_fdr <- p.adjust(p = cor_stats_pro_tab$genepair_cor_p_value, method = "fdr")


# test correlation with phopho level of other genes -----------------------
row_indexs <- which(rowSums(!is.na(phospho_tab[,geneA_exp_nonna_ids])) >= least_no_datapoints)

sample_sizes <- vector(mode = "numeric", length = length(row_indexs))
p_values <- vector(mode = "numeric", length = length(row_indexs))
rhos <- vector(mode = "numeric", length = length(row_indexs))

for (i in 1:length(row_indexs)) {
  geneB_exp_values <- phospho_tab[row_indexs[i], tumor_bulk_aliquot_ids]
  
  nonna_ids <- tumor_bulk_aliquot_ids[!is.na(geneA_exp_values) & !is.na(geneB_exp_values)]
  
  cor_stats_tmp <- cor.test(x = unlist(geneA_exp_values[,nonna_ids]), y = unlist(geneB_exp_values[,nonna_ids]), method = "spearman")
  p_values[i] <- cor_stats_tmp$p.value
  rhos[i] <- cor_stats_tmp$estimate
  sample_sizes[i] <- length(nonna_ids)
  print(i)
}

cor_stats_phospho_tab <- data.frame(geneB = phospho_tab[row_indexs, "Gene"],
                                    geneB_exp_type = "PHO", 
                                    geneB_exp_site = str_split_fixed(string = phospho_tab[row_indexs, "Index"], pattern = "_", n = 7)[,7],
                                    genepair_cor_p_value = p_values,
                                    genepair_cor_rho = rhos,
                                    genepair_cor_sample_size = sample_sizes,
                                    geneA = geneA, 
                                    geneA_exp_type = geneA_exp_type)
cor_stats_phospho_tab$genepair_cor_fdr <- p.adjust(p = cor_stats_phospho_tab$genepair_cor_p_value, method = "fdr")

# merge all the result table ----------------------------------------------
cor_stats_sup_tab <- rbind(cor_stats_rna_tab, cor_stats_pro_tab, cor_stats_phospho_tab)
write.table(x = cor_stats_sup_tab, file = paste0(dir_out, "Spearman_Cor_Stats_with.", geneA, ".", geneA_exp_type, ".tsv"), quote = F, sep = "\t", row.names = F)




