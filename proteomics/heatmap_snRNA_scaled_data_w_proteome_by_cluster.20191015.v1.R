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
sn_ave_scaled_exp_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/proteomics/test_cor_sn_cluster_average_w_proteome/20191015.v1/Average_Scaled_Data_for_Integrated.20191015.v1.Run.20191015.v1.tsv", data.table = F)
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


# input spearman correlation results --------------------------------------
spearman.test.result.df.sup <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/proteomics/test_cor_sn_cluster_average_w_proteome/20191015.v1/Spearman_Correlation_btw_Average_Scaled_snRNA_data_w_bulk_Protein_for_Integrated.20191015.v1.Run.20191015.v1.tsv", data.table = F)

# merge data  -------------------------------------------------------------
for (gene_tmp in unique(spearman.test.result.df.sup$gene_symbol[spearman.test.result.df.sup$p.value < 0.05 & spearman.test.result.df.sup$rho > 0])) {
  pro.m <- melt(pro_mat[gene_tmp, meta_tab_in_use$Specimen.ID.snRNA])
  pro.m$aliquot <- rownames(pro.m)
  pro.m <- pro.m %>%
    dplyr::rename(exp_value = value) %>%
    mutate(column_text = "Protein") %>%
    mutate(cluster = NA)
  
  snRNA.m <- melt(sn_ave_scaled_exp_tab[gene_tmp,])
  snRNA.m$cluster_aliquot <- rownames(snRNA.m)
  snRNA.m <- snRNA.m %>%
    mutate(cluster = str_split_fixed(string = cluster_aliquot, pattern = "_", n = 2)[,1]) %>%
    mutate(aliquot = str_split_fixed(string = cluster_aliquot, pattern = "_", n = 2)[,2]) %>%
    filter(aliquot %in% meta_tab_in_use$Specimen.ID.snRNA) %>%
    mutate(column_text = paste0("snRNA_C", cluster)) %>%
    dplyr::rename(exp_value = value) %>%
    select(exp_value, aliquot, column_text, cluster)
  
  # plot scatterplots -------------------------------------------------------
  tab2p <- merge(pro.m, snRNA.m, by = c("aliquot"), suffixes = c(".pro", ".snRNA"))
  spearman.test.result.tmp <- spearman.test.result.df.sup %>%
    filter(gene_symbol == gene_tmp)
  
  tab2p <- merge(tab2p, spearman.test.result.tmp, by.x = c("cluster.snRNA"), by.y = c("cluster"), all.x = T)
  tab2p_line_highlight <- tab2p %>%
    filter(p.value < 0.05 & rho > 0)
  
  tab2p_line_vague <- tab2p %>%
    filter(!(p.value < 0.05 & rho > 0))
  
  p <- ggplot()
  p <- p + geom_point(data = tab2p, mapping = aes(x = exp_value.pro , y = exp_value.snRNA, color = column_text.snRNA), size = 2)
  p <- p + geom_line(data = tab2p_line_vague, mapping = aes(x = exp_value.pro , y = exp_value.snRNA, color = column_text.snRNA), alpha = 0.2)
  p <- p + geom_line(data = tab2p_line_highlight, mapping = aes(x = exp_value.pro , y = exp_value.snRNA, color = column_text.snRNA))
  p <- p + theme_bw()
  file2write <- paste(dir_out, "Scatterplot_Protein_vs_snRNA_", gene_tmp, ".", run_id, ".png", sep="")
  png(file = file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
}



# plot heatmap ------------------------------------------------------------
for (gene_tmp in c("CDH13", "VHL")) {
  pro.m <- melt(pro_mat[gene_tmp, meta_tab_in_use$Specimen.ID.snRNA])
  pro.m$aliquot <- rownames(pro.m)
  pro.m <- pro.m %>%
    dplyr::rename(exp_value = value) %>%
    mutate(column_text = "Protein") %>%
    mutate(cluster = NA)
  
  snRNA.m <- melt(sn_ave_scaled_exp_tab[gene_tmp,])
  snRNA.m$cluster_aliquot <- rownames(snRNA.m)
  snRNA.m <- snRNA.m %>%
    mutate(cluster = str_split_fixed(string = cluster_aliquot, pattern = "_", n = 2)[,1]) %>%
    mutate(aliquot = str_split_fixed(string = cluster_aliquot, pattern = "_", n = 2)[,2]) %>%
    filter(aliquot %in% meta_tab_in_use$Specimen.ID.snRNA) %>%
    mutate(column_text = paste0("snRNA_C", cluster)) %>%
    dplyr::rename(exp_value = value) %>%
    select(exp_value, aliquot, column_text, cluster)
  
  
  tab2p <- rbind(pro.m, snRNA.m)
  p <- ggplot()
  p <- p + geom_tile(data = tab2p, mapping = aes(x = column_text , y = aliquot, fill = exp_value))
  p <- p + scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab" )
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5))
  file2write <- paste(dir_out, "Heatmap_Protein_vs_snRNA_", gene_tmp, ".", run_id, ".png", sep="")
  png(file = file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
}


