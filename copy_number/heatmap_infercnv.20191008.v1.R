# Yige Wu @WashU Sep 2019
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
dir_infercnv_out <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/inferCNV/outputs/integration.20190927.v1/"

# set samples to process ----------------------------------------------------------
sample_ids <- c("CPT0075130004", "CPT0001260013", "CPT0075140002", "CPT0086350004")


# set genes 2 plot --------------------------------------------------------
genes_loss <- c("VHL", "PBRM1", "BAP1", "SETD2",
                "HIF1A",
                "CDKN2A", "PTEN", 
                "NEGR1",
                "QKI",
                "CADM2", 
                "PTPRD", "NRXN3")
genes_gain <- c("PRKCI", 
                "MECOM",
                "MDM4",
                "MYC",
                "JAK2")
ccrcc_cna_genes_df <- data.frame(gene_symbol = c(genes_loss, genes_gain),
                                 gene_cna_type = c(rep("Loss", length(genes_loss)), rep("Gain", length(genes_gain))))
ccrcc_cna_genes_df

# input infercnv results ------------------------------------------------
sup_tab_infercnv <- NULL
for (sample_id_tmp in sample_ids) {
  dir_infercnv_tmp <- paste0(dir_infercnv_out, sample_id_tmp, "/")
  path_hmmi6_tumor_file_tmp <- paste0(dir_infercnv_tmp, "infercnv.14_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.observations.txt")
  hmmi6_tumor_tab_tmp <- fread(input = path_hmmi6_tumor_file_tmp, data.table = F)
  number_tumor_cells <- ncol(hmmi6_tumor_tab_tmp) - 1
  number_tumor_cells
  
  path_hmmi6_normal_file_tmp <- paste0(dir_infercnv_tmp, "infercnv.14_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.references.txt")
  hmmi6_normal_tab_tmp <- fread(input = path_hmmi6_normal_file_tmp, data.table = F)
  number_normal_cells <- ncol(hmmi6_normal_tab_tmp) - 1
  number_normal_cells
  
  hmmi6_tab_tmp <- hmmi6_tumor_tab_tmp %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% ccrcc_cna_genes_df$gene_symbol)
  
  hmmi6_tab_m_tmp <- melt(hmmi6_tab_tmp, id.vars = c("gene_symbol"))
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
    mutate(aliquot = sample_id_tmp)
  
  sup_tab_infercnv <- rbind(hmmi6_tab_m_tmp, sup_tab_infercnv)
}
sup_tab_infercnv <- sup_tab_infercnv %>%
  mutate(cna_text = paste0(gene_symbol, "_", gene_cna_type)) %>%
  mutate(cna_cat = map_infercnv_state2category(copy_state = cna_state))


# load meta data ----------------------------------------------------------
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)
sup_tab_infercnv <- merge(sup_tab_infercnv, meta_tab, by.x = c("aliquot"), by.y = c("Specimen.ID.snRNA"), all.x = T)
sup_tab_infercnv$Case.ID[sup_tab_infercnv$aliquot == "CPT0075140002"] <- "C3N-01200"

# plot heatmap ------------------------------------------------------------
tab2p <- sup_tab_infercnv

p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = cna_text, y = aliquot, size = perc_cna_in_tumor_cell, fill = cna_cat), shape = 21)
p <- p + scale_fill_manual(values = cna_state_colors)
p <- p + facet_grid(Case.ID~gene_cna_type, scales = "free", space = "free", shrink = T)
p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p

png(file = paste0(dir_out, "Percent_Tumor_with_CNA", ".", run_id, ".png"), width = 1000, height = 600, res = 150)
print(p)
dev.off()
