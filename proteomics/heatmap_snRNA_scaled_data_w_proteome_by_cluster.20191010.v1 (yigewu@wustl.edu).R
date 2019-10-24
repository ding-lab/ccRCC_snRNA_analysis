# Yige Wu @WashU Oct 2019
## for plotting the scaled snRNA expression with bulk protein to show correlation across cell clusters

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# specify genes to plot ---------------------------------------------------
genes2plot <- c("KRT7")

# input scaled snRNA by cluster -------------------------------------------
sn_ave_scaled_exp_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/proteomics/test_cor_sn_cluster_average_w_proteome/20191010.v2/Average_Scaled_Data_for_Integrated_5FACS.20190927.v1.Run.20191010.v1.tsv", data.table = F)
## format to matrix
rownames(sn_ave_scaled_exp_tab) <- sn_ave_scaled_exp_tab$V1
sn_ave_scaled_exp_tab$V1 <- NULL
sn_ave_scaled_exp_tab <- as.matrix(sn_ave_scaled_exp_tab)
aliquot_ids <- colnames(sn_ave_scaled_exp_tab)
aliquot_ids <- unique(str_split_fixed(string = aliquot_ids, pattern = "_", n = 2)[,2])
aliquot_ids

# input protein data ------------------------------------------------------
pro_tab <- loadParseProteomicsData(expression_type  = "PRO", sample_type = "tumor")
pro_mat <- pro_tab
rownames(pro_mat) <- pro_tab$Gene
pro_mat$Gene <- NULL
pro_mat <- as.matrix(pro_mat)


# reformat protein data colnames to aliquot ids ---------------------------
## input meta data
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)

meta_tab_in_use <- meta_tab %>%
  dplyr::filter(meta_tab$Specimen.ID.snRNA %in% aliquot_ids)
rownames(meta_tab_in_use) <- meta_tab_in_use$Case.ID

pro_mat <- pro_mat[,meta_tab_in_use$Case.ID]
pro_mat
colnames(pro_mat) <- meta_tab_in_use$Specimen.ID.snRNA

# merge data  -------------------------------------------------------------
gene_tmp <- "KRT7"

pro.m <- melt(pro_mat[gene_tmp, meta_tab_in_use$Specimen.ID.snRNA])
pro.m$aliquot <- rownames(pro.m)
pro.m <- pro.m %>%
  dplyr::rename(exp_value = value) %>%
  mutate(col_text = "Bulk_Protein")

snRNA.m <- melt(sn_ave_scaled_exp_tab[gene_tmp,])
snRNA.m$cluster_aliquot <- rownames(snRNA.m)
snRNA.m <- snRNA.m %>%
  mutate(cluster = str_split_fixed(string = cluster_aliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(aliquot = str_split_fixed(string = cluster_aliquot, pattern = "_", n = 2)[,2]) %>%
  dplyr::rename(snRNA_scaled_value = value)


tab2p <- merge(pro.m, snRNA.m, by = c("aliquot"), all.y = T)

p <- ggplot()
p <- p + geom_tile(data = tab2p, mapping = aes())

tab2p


