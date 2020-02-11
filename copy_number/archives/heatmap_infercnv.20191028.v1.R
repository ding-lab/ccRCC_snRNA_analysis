# Yige Wu @WashU Oct 2019
## for plotting the genomics landscape of integrated single cell data

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# set infercnv output directory -------------------------------------------
dir_infercnv_output <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/"
integration_id <- "20191021.v1"

# set samples to process ----------------------------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# input infercnv results ------------------------------------------------
sup_tab_infercnv <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  tumor_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  number_tumor_cells <- ncol(tumor_cnv_state_mat) - 1
  number_tumor_cells
  
  ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  number_normal_cells <- ncol(ref_cnv_state_mat) - 1
  number_normal_cells
  
  infercnv_obj <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  
  
  tumor_cnv_state_mat <- tumor_cnv_state_mat %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% ccrcc_cna_genes_df$gene_symbol)
  
  ref_cnv_state_mat <- ref_cnv_state_mat %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% ccrcc_cna_genes_df$gene_symbol)
  
  tumor_cnv_state_mat.m <- melt(tumor_cnv_state_mat, id.vars = c("gene_symbol"))
  ref_cnv_state_mat.m <- melt(ref_cnv_state_mat, id.vars = c("gene_symbol"))
  
  cnv_state_mat.m <- rbind(tumor_cnv_state_mat.m, ref_cnv_state_mat.m)
  cnv_state_mat.m <- cnv_state_mat.m %>%
    rename(cna_state = value) %>%
    filter(cna_state != 1) %>%
    select(gene_symbol, cna_state) %>%
    table() %>%
    as.data.frame()
  cnv_state_mat.m <- merge(cnv_state_mat.m, ccrcc_cna_genes_df, by = c("gene_symbol"))
  cnv_state_mat.m <- cnv_state_mat.m %>%
    filter((gene_cna_type == "Loss" & cna_state %in% c(0.5, 0)) | (gene_cna_type == "Gain" & cna_state %in% c(1.5, 2))) %>%
    mutate(num_tumor_cell = number_tumor_cells) %>%
    mutate(num_all_cell = (number_tumor_cells + number_normal_cells)) %>%
    mutate(perc_cna_in_tumor_cell = (Freq/number_tumor_cells)) %>%
    mutate(perc_cna_in_all_cell = (Freq/(number_tumor_cells + number_normal_cells))) %>%
    mutate(aliquot = snRNA_aliquot_id_tmp)
  
  sup_tab_infercnv <- rbind(cnv_state_mat.m, sup_tab_infercnv)
}
sup_tab_infercnv <- sup_tab_infercnv %>%
  mutate(cna_text = paste0(gene_symbol, "_", gene_cna_type)) %>%
  mutate(sn_cna_cat = map_infercnv_state2category(copy_state = cna_state))

# load meta data ----------------------------------------------------------
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)
sup_tab <- merge(sup_tab_infercnv, meta_tab, by.x = c("aliquot"), by.y = c("Specimen.ID.snRNA"), all.x = T)
sup_tab$Case.ID[sup_tab$aliquot == "CPT0075140002"] <- "C3N-01200"
sup_tab$Case.ID[sup_tab$aliquot == "CPT0025890002"] <- "C3N-00733"


# add chromosomal region --------------------------------------------------
sup_tab$chr_region <- mapvalues(x = sup_tab$gene_symbol, from = ccrcc_cna_genes_df$gene_symbol, to = as.vector(ccrcc_cna_genes_df$chr_region))
sup_tab$chr_region <- factor(sup_tab$chr_region, levels = c("3p", "5q", "14q",
                                                            "9p21", 
                                                            "10q23",
                                                            "1p31",
                                                            "6q24",
                                                            "3p12",
                                                            "9p23",
                                                            "14q24",
                                                            "3p26",
                                                            "1q32",
                                                            "8q24",
                                                            "9q24"))

# plot heatmap showing sn-based CNA Frequency------------------------------------------------------------
tab2p <- sup_tab
tab2p$perc_cna_in_all_cell[tab2p$perc_cna_in_all_cell == 0] <- NA

p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = gene_symbol, y = aliquot, size = perc_cna_in_all_cell, fill = sn_cna_cat), shape = 21)
p <- p + scale_fill_manual(values = copy_number_colors)
p <- p + facet_grid(Case.ID~chr_region + gene_cna_type, scales = "free", space = "free", shrink = T)
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p <- p + theme(panel.spacing = unit(0, "lines"))
p
png(file = paste0(dir_out, "TCGA_CNA_Genes_snCNA.withlegend", ".", run_id, ".png"), width = 1200, height = 800, res = 150)
print(p)
dev.off()

# write out CNA frequency table -------------------------------------------
write.table(x = tab2p, file = paste0(dir_out, "TCGA_CNA_Genes_snCNA_Frequency.tsv"), quote = F, sep = "\t", row.names = F)
