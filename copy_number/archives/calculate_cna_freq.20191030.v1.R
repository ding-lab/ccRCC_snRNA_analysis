# Yige Wu @WashU Oct 2019
## for plotting the genomics landscape of integrated single cell data

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
packages = c(
  "GWASTools"
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
  }
  library(package = pkg_name_tmp, character.only = T)
}

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


# set genes to look for CNA -----------------------------------------------
cna_genes <- unique(c(as.vector(ccrcc_cna_genes_df$gene_symbol), SMGs[["CCRCC"]]))
cna_genes


# calculate cna frequency by aliquot --------------------------------------
gene_cna_state_tab <- NULL
number_cells_tab <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  
  # calculate cna frequency by gene --------------------------------------
  tumor_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  number_tumor_cells <- ncol(tumor_cnv_state_mat) - 1
  number_tumor_cells
  
  ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  number_normal_cells <- ncol(ref_cnv_state_mat) - 1
  number_normal_cells
  
  number_cells_tab <- rbind(number_cells_tab, data.frame(snRNA_aliquot_id = snRNA_aliquot_id_tmp, number_tumor_cells = number_tumor_cells, number_normal_cells = number_normal_cells))
  
  tumor_cnv_state_mat <- tumor_cnv_state_mat %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% cna_genes)
  
  ref_cnv_state_mat <- ref_cnv_state_mat %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% cna_genes)
  
  tumor_cnv_state_mat.m <- melt(tumor_cnv_state_mat, id.vars = c("gene_symbol"))
  ref_cnv_state_mat.m <- melt(ref_cnv_state_mat, id.vars = c("gene_symbol"))
  
  cnv_state_mat.m <- rbind(tumor_cnv_state_mat.m, ref_cnv_state_mat.m)
  cnv_state_mat.m <- cnv_state_mat.m %>%
    rename(cna_state = value) %>%
    select(gene_symbol, cna_state) %>%
    table() %>%
    as.data.frame()
  cnv_state_mat.m <- cnv_state_mat.m %>%
    mutate(num_tumor_cell = number_tumor_cells) %>%
    mutate(num_all_cell = (number_tumor_cells + number_normal_cells)) %>%
    mutate(perc_cna_in_tumor_cell = (Freq/number_tumor_cells)) %>%
    mutate(perc_cna_in_all_cell = (Freq/(number_tumor_cells + number_normal_cells))) %>%
    mutate(snRNA_aliquot_id = snRNA_aliquot_id_tmp)
  
  gene_cna_state_tab <- rbind(cnv_state_mat.m, gene_cna_state_tab)
}

# write out CNA frequency by gene table -------------------------------------------
write.table(x = gene_cna_state_tab, file = paste0(dir_out, "TCGA_CNA_Genes_snCNA_Frequency.tsv"), quote = F, sep = "\t", row.names = F)

# calculate cna frequency by chr arm --------------------------------------
## input centromere positions
data(centromeres.hg38)

## create super table to store the cna states
chr_arm_cna_state_tab <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  cnv_region_tab <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"), data.table = F)
  
  ## input cell group 2 cell barcode table and estimate the number of cells in each cell group
  cell_groupings_tab <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.cell_groupings"), data.table = F)
  cell_groupings_counts_tab <- table(cell_groupings_tab$cell_group_name)
  cell_groupings_counts_tab <- as.data.frame(cell_groupings_counts_tab)
  
  ## estimate 3p loss frequency
  chr3p_loss_tab <- cnv_region_tab %>%
    filter(chr == "chr3") %>%
    mutate(is_3p = (start < centromeres.hg38$left.base[centromeres.hg38$chrom == 3])) %>%
    filter(is_3p) %>%
    filter(state < 1) %>%
    select(cell_group_name) %>%
    unique() %>%
    mutate(cell_group_cell_count = as.numeric(mapvalues(x = cell_group_name, from = cell_groupings_counts_tab$Var1, to = as.vector(cell_groupings_counts_tab$Freq))))
  
  chr3p_loss_count <- chr3p_loss_tab %>%
    summarize(Freq = sum(cell_group_cell_count, na.rm = T)) %>%
    mutate(chr_region = "3p") %>%
    mutate(gene_cna_type = "Loss")
  
  ## estimate 5q gain frequency
  chr5q_gain_tab <- cnv_region_tab %>%
    filter(chr == "chr5") %>%
    mutate(is_5q = (end > centromeres.hg38$right.base[centromeres.hg38$chrom == 5])) %>%
    filter(is_5q) %>%
    filter(state > 1) %>%
    select(cell_group_name) %>%
    unique() %>%
    mutate(cell_group_cell_count = as.numeric(mapvalues(x = cell_group_name, from = cell_groupings_counts_tab$Var1, to = as.vector(cell_groupings_counts_tab$Freq))))
  chr5q_gain_count <- chr5q_gain_tab %>%
    summarize(Freq = sum(cell_group_cell_count, na.rm = T)) %>%
    mutate(chr_region = "5q") %>%
    mutate(gene_cna_type = "Gain")
  
  ## estimate 14q loss frequency
  chr14q_loss_tab <- cnv_region_tab %>%
    filter(chr == "chr14") %>%
    mutate(is_14q = (end > centromeres.hg38$right.base[centromeres.hg38$chrom == 14])) %>%
    filter(is_14q) %>%
    filter(state < 1) %>%
    select(cell_group_name) %>%
    unique() %>%
    mutate(cell_group_cell_count = as.numeric(mapvalues(x = cell_group_name, from = cell_groupings_counts_tab$Var1, to = as.vector(cell_groupings_counts_tab$Freq))))
  
  chr14q_loss_count <- chr14q_loss_tab %>%
    summarize(Freq = sum(cell_group_cell_count, na.rm = T)) %>%
    mutate(chr_region = "14q") %>%
    mutate(gene_cna_type = "Loss")
  
  chr_arm_cna_state_tmp <- rbind(chr3p_loss_count, chr5q_gain_count, chr14q_loss_count)
  chr_arm_cna_state_tmp$snRNA_aliquot_id <- snRNA_aliquot_id_tmp
  chr_arm_cna_state_tab <- rbind(chr_arm_cna_state_tab, chr_arm_cna_state_tmp)
}
chr_arm_cna_state_tab <- merge(chr_arm_cna_state_tab, number_cells_tab, by = c("snRNA_aliquot_id"), all.x = T)
chr_arm_cna_state_tab <- chr_arm_cna_state_tab %>%
  mutate(num_all_cell = (number_tumor_cells + number_normal_cells)) %>%
  mutate(perc_cna_in_tumor_cell = (Freq/number_tumor_cells)) %>%
  mutate(perc_cna_in_all_cell = (Freq/(number_tumor_cells + number_normal_cells)))
  
write.table(x = chr_arm_cna_state_tab, file = paste0(dir_out, "TCGA_CNA_Chr_Arm_snCNA_Frequency.tsv"), quote = F, sep = "\t", row.names = F)

