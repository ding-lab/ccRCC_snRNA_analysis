# Yige Wu @WashU Feb 2020
## for calculating the fraction of tumor cells with cnv in different frequently altered genes

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(dplyr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## set infercnv output directory
dir_infercnv_output <- "./Resources/snRNA_Processed_Data/InferCNV/outputs/"
infercnv_run_id <- "Individual.20200305.v1"
dir_infercnv_run <- paste0(dir_infercnv_output, infercnv_run_id, "/")
## get aliquots to process
aliquots2process <- list.files(dir_infercnv_run)
aliquots2process <- aliquots2process[grepl(pattern = "CPT", x = aliquots2process)]
## input barcode to cell type info
barcode2cluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20201130.v1/Barcode2TumorSubclusterId.20201130.v1.tsv", data.table = F)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")

# input infercnv results ------------------------------------------------
cnv_state_count_aliquots <- NULL
for (aliquot_tmp in aliquots2process) {
  ## input infercnv cnv state outputs
  cnv_state_obs_mat <- fread(input = paste0(dir_infercnv_run, aliquot_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  cnv_state_ref_mat <- fread(input = paste0(dir_infercnv_run, aliquot_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  
  ## dplyr::filter cnv results to only selected genes
  cnv_state_obs_mat <- cnv_state_obs_mat %>%
    dplyr::rename(gene_symbol = V1) %>%
    dplyr::filter(gene_symbol %in% knowncnvgenes_df$Gene_Symbol)
  cnv_state_ref_mat <- cnv_state_ref_mat %>%
    dplyr::rename(gene_symbol = V1) %>%
    dplyr::filter(gene_symbol %in% knowncnvgenes_df$Gene_Symbol)
  
  ## melt wide data frame to long data frame
  cnv_state_obs_mat.m <- melt(cnv_state_obs_mat, id.vars = c("gene_symbol"))
  cnv_state_ref_mat.m <- melt(cnv_state_ref_mat, id.vars = c("gene_symbol"))
  
  ## combine the cnv results
  cnv_state_mat.m <- rbind(cnv_state_obs_mat.m, cnv_state_ref_mat.m)
  
  ## get the barcodes for this aliquot
  aliquot_barcode2cluster_df <- barcode2cluster_df %>%
    dplyr::filter(orig.ident == aliquot_tmp) %>%
    dplyr::mutate(Id_TumorManualCluster = Cluster_Name) %>%
    dplyr::mutate(individual_barcode = barcode)
  
  ## dplyr::rename columns and dplyr::filter down to only maligant nephron epithlium cells
  cnv_state_mat.m <- cnv_state_mat.m %>%
    dplyr::rename(cna_state = value) %>%
    dplyr::rename(barcode = variable) %>%
    dplyr::filter(barcode %in% aliquot_barcode2cluster_df$individual_barcode)
  ## map barcode to tumor subcluster
  colnames(barcode2cluster_df)
  cnv_state_mat.m$tumor_subcluster <- mapvalues(x = cnv_state_mat.m$barcode, from = aliquot_barcode2cluster_df$individual_barcode, to = aliquot_barcode2cluster_df$Id_TumorManualCluster)
  
  ## count number of cells with different cnv state per gene
  cnv_state_count <- cnv_state_mat.m %>%
    dplyr::select(gene_symbol, tumor_subcluster, cna_state) %>%
    table() %>%
    as.data.frame() %>%
    dplyr::filter(Freq > 0) %>%
    dplyr::mutate(aliquot = aliquot_tmp)
  
  ## summarize the detected values for each gene
  cnv_nonna_count <- cnv_state_count %>%
    group_by(gene_symbol, tumor_subcluster) %>%
    summarize(num_cells_nonna = sum(Freq))
  ## merge counts
  cnv_state_count <- merge(cnv_state_count, cnv_nonna_count, by = c("gene_symbol", "tumor_subcluster"), all.x = T)
  
  ## combine with super table
  cnv_state_count_aliquots <- rbind(cnv_state_count, cnv_state_count_aliquots)
}

# write out CNA frequency by 6 state by cluster -------------------------------------------
## calculate fraction
tmp <- cnv_state_count_aliquots
tmp$Fraction <- tmp$Freq/(tmp$num_cells_nonna)
tmp$cna_state <- as.vector(tmp$cna_state)
cnv_6state_count_aliquots <- tmp
## write
write.table(x = cnv_6state_count_aliquots, file = paste0(dir_out, "fraction_of_tumorcells_with_cnv_by_gene_by_6state.per_manualsubcluster.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

# write out CNA frequency by 3 state by cluster -------------------------------------------
## calculate fraction
tmp <- cnv_state_count_aliquots
tmp$Fraction <- tmp$Freq/(tmp$num_cells_nonna)
tmp$cna_state <- as.numeric(as.vector(tmp$cna_state))
## annotate cells with expected cnv state
cnv_3state_count_aliquots <- tmp %>%
  dplyr::mutate(cna_3state = ifelse(cna_state %in% c(1.5, 2, 3), "Gain", ifelse(cna_state %in% c(0, 0.5), "Loss", "Neutral"))) %>%
  group_by(aliquot, tumor_subcluster, gene_symbol, cna_3state) %>%
  summarise(Fraction = sum(Fraction, na.rm = T))
## write
write.table(x = cnv_3state_count_aliquots, file = paste0(dir_out, "fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)



