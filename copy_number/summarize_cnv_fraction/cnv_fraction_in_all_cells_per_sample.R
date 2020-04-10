# Yige Wu @WashU Feb 2020
## for plotting the fraction of cells with CNV per sample in case per cluster CNV distribution is too confusing

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
## load meta data
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## set infercnv output directory
dir_infercnv_output <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/"
infercnv_run_id <- "Individual.20200207.v1"
dir_infercnv_run <- paste0(dir_infercnv_output, infercnv_run_id, "/")
## get aliquots to process
aliquots2process <- list.files(dir_infercnv_run)

# input infercnv results ------------------------------------------------
cnv_state_count_aliquots <- NULL
for (aliquot_tmp in aliquots2process) {
  ## input infercnv cnv state outputs
  cnv_state_obs_mat <- fread(input = paste0(dir_infercnv_run, aliquot_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  cnv_state_ref_mat <- fread(input = paste0(dir_infercnv_run, aliquot_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  
  ## get the number of cells
  number_tumor_cells <- ncol(cnv_state_obs_mat) - 1
  number_tumor_cells
  number_normal_cells <- ncol(cnv_state_ref_mat) - 1
  number_normal_cells
  
  ## filter cnv results to only selected genes
  cnv_state_obs_mat <- cnv_state_obs_mat %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% ccrcc_cna_genes_df$gene_symbol)
  cnv_state_ref_mat <- cnv_state_ref_mat %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% ccrcc_cna_genes_df$gene_symbol)
  
  ## melt wide data frame to long data frame
  cnv_state_obs_mat.m <- melt(cnv_state_obs_mat, id.vars = c("gene_symbol"))
  cnv_state_ref_mat.m <- melt(cnv_state_ref_mat, id.vars = c("gene_symbol"))
  
  ## combine the cnv results
  cnv_state_mat.m <- rbind(cnv_state_obs_mat.m, cnv_state_ref_mat.m)
  
  ## count number of cells with different cnv state per gene
  cnv_state_count <- cnv_state_mat.m %>%
    rename(cna_state = value) %>%
    select(gene_symbol, cna_state) %>%
    table() %>%
    as.data.frame() %>%
    filter(Freq > 0)
  
  ## annotate each gene with expected cnv state and chromosome region
  cnv_state_count <- merge(cnv_state_count, ccrcc_cna_genes_df, by = c("gene_symbol"))
  
  ## annotate cells with expected cnv state
  cnv_state_count <- cnv_state_count %>%
    mutate(expected_cna = (gene_cna_type == "Loss" & cna_state %in% c(0.5, 0)) | (gene_cna_type == "Gain" & cna_state %in% c(1.5, 2))) %>%
    mutate(aliquot = aliquot_tmp)
  
  ## summarize the detected values for each gene
  cnv_nonna_count <- cnv_state_count %>%
    group_by(gene_symbol) %>%
    summarize(num_cells_nonna = sum(Freq))
  
  ## merge counts
  cnv_state_count <- merge(cnv_state_count, cnv_nonna_count, by = c("gene_symbol"), all.x = T)
  
  ## calculate fraction
  cnv_state_count$Fraction <- cnv_state_count$Freq/(cnv_state_count$num_cells_nonna)
  
  ## combine with super table
  cnv_state_count_aliquots <- rbind(cnv_state_count, cnv_state_count_aliquots)
}
cnv_state_count_aliquots <- cnv_state_count_aliquots %>%
  mutate(cna_text = paste0(gene_symbol, "_", gene_cna_type)) %>%
  mutate(sn_cna_cat = map_infercnv_state2category(copy_state = cna_state))

# write out CNA frequency table -------------------------------------------
write.table(x = cnv_state_count_aliquots, file = paste0(dir_out, "TCGA_CNA_Genes_snCNA_Frequency.tsv"), quote = F, sep = "\t", row.names = F)
