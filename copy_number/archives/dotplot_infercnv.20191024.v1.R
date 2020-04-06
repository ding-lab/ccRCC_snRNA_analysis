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
dir_infercnv_out <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/integration.20191021.v1/"

# set samples to process ----------------------------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# input infercnv results ------------------------------------------------
sup_tab_infercnv <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  dir_infercnv_tmp <- paste0(dir_infercnv_out, snRNA_aliquot_id_tmp, "/")
  
  path_hmmi6_tumor_file_tmp <- paste0(dir_infercnv_tmp, "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt")
  hmmi6_tumor_tab_tmp <- fread(input = path_hmmi6_tumor_file_tmp, data.table = F)
  number_tumor_cells <- ncol(hmmi6_tumor_tab_tmp) - 1
  number_tumor_cells
  
  path_hmmi6_normal_file_tmp <- paste0(dir_infercnv_tmp, "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt")
  hmmi6_normal_tab_tmp <- fread(input = path_hmmi6_normal_file_tmp, data.table = F)
  number_normal_cells <- ncol(hmmi6_normal_tab_tmp) - 1
  number_normal_cells
  
  hmmi6_tumor_tab_tmp <- hmmi6_tumor_tab_tmp %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% ccrcc_cna_genes_df$gene_symbol)
  
  hmmi6_normal_tab_tmp <- hmmi6_normal_tab_tmp %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% ccrcc_cna_genes_df$gene_symbol)
  
  hmmi6_tumor_tab_m_tmp <- melt(hmmi6_tumor_tab_tmp, id.vars = c("gene_symbol"))
  hmmi6_normal_tab_m_tmp <- melt(hmmi6_normal_tab_tmp, id.vars = c("gene_symbol"))
  
  hmmi6_tab_m_tmp <- rbind(hmmi6_tumor_tab_m_tmp, hmmi6_normal_tab_m_tmp)
  hmmi6_tab_m_tmp <- hmmi6_tab_m_tmp %>%
    rename(cna_state = value) %>%
    filter(cna_state != 1) %>%
    select(gene_symbol, cna_state) %>%
    table() %>%
    as.data.frame()
  hmmi6_tab_m_tmp <- merge(hmmi6_tab_m_tmp, ccrcc_cna_genes_df, by = c("gene_symbol"))
  hmmi6_tab_m_tmp <- hmmi6_tab_m_tmp %>%
    filter((gene_cna_type == "Loss" & cna_state %in% c(0.5, 0)) | (gene_cna_type == "Gain" & cna_state %in% c(1.5, 2))) %>%
    mutate(num_tumor_cell = number_tumor_cells) %>%
    mutate(num_all_cell = (number_tumor_cells + number_normal_cells)) %>%
    mutate(perc_cna_in_tumor_cell = (Freq/number_tumor_cells)) %>%
    mutate(perc_cna_in_all_cell = (Freq/(number_tumor_cells + number_normal_cells))) %>%
    mutate(aliquot = snRNA_aliquot_id_tmp)
  
  sup_tab_infercnv <- rbind(hmmi6_tab_m_tmp, sup_tab_infercnv)
}
sup_tab_infercnv <- sup_tab_infercnv %>%
  mutate(cna_text = paste0(gene_symbol, "_", gene_cna_type)) %>%
  mutate(sn_cna_cat = map_infercnv_state2category(copy_state = cna_state))


# load meta data ----------------------------------------------------------
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)
sup_tab <- merge(sup_tab_infercnv, meta_tab, by.x = c("aliquot"), by.y = c("Specimen.ID.snRNA"), all.x = T)
sup_tab$Case.ID[sup_tab$aliquot == "CPT0075140002"] <- "C3N-01200"

# plot heatmap showing sn-based CNA and somatic mutation and VAF------------------------------------------------------------
tab2p <- sup_tab
tab2p$perc_cna_in_all_cell[tab2p$perc_cna_in_all_cell == 0] <- NA
tab2p <- tab2p %>%
  mutate(vaf_summary_text = paste0(signif(x = perc_cna_in_all_cell*100, digits = 2), "% cells w. CNA\n", 
                               signif(x = vaf*100, digits = 2), "% WES reads w. variant")) %>%
  mutate(ggrepel_text = ifelse(is.na(vaf), NA, vaf_summary_text))

p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = gene_symbol, y = aliquot, size = perc_cna_in_all_cell, fill = sn_cna_cat), shape = 21)
p <- p + scale_fill_manual(values = copy_number_colors)
p <- p + geom_point(data = tab2p, mapping = aes(x = gene_symbol, y = aliquot, shape = somatic_variant_class, size = vaf))
p <- p + scale_shape_manual(values = c("Frame_Shift_Del" = 7, "Frame_Shift_Ins" = 7, "Nonsense_Mutation" = 7, "Missense_Mutation" = 12))
p <- p + facet_grid(Case.ID~gene_cna_type, scales = "free", space = "free", shrink = T)
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p
png(file = paste0(dir_out, "TCGA_CNA_Genes_snCNA_bulkMutation.withlegend", ".", run_id, ".png"), width = 1200, height = 800, res = 150)
print(p)
dev.off()

