# Yige Wu @WashU Feb 2020
## for calculating the fraction of tumor cells with cnv in different frequently altered genes

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
## set infercnv output directory
dir_infercnv_output <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/"
infercnv_run_id <- "Individual.20200207.v1"
dir_infercnv_run <- paste0(dir_infercnv_output, infercnv_run_id, "/")
## get aliquots to process
aliquots2process <- list.files(dir_infercnv_run)
## input barcode to cell type info
barcode2cluster_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/fetch_data_malignant_nephron_epithelium_cells_reclustered/20200225.v1/Cluster_UMAP_Data.Malignant_Nephron_Epithelium.20200225.v1.tsv", data.table = F)

# input infercnv results ------------------------------------------------
cnv_state_count_aliquots <- NULL
for (aliquot_tmp in aliquots2process) {
  ## input infercnv cnv state outputs
  cnv_state_obs_mat <- fread(input = paste0(dir_infercnv_run, aliquot_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  cnv_state_ref_mat <- fread(input = paste0(dir_infercnv_run, aliquot_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  
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
  
  ## get the barcodes for this aliquot
  aliquot_barcode2cluster_df <- barcode2cluster_df %>%
    filter(orig.ident == aliquot_tmp)
  
  ## rename columns and filter down to only maligant nephron epithlium cells
  cnv_state_mat.m <- cnv_state_mat.m %>%
    rename(cna_state = value) %>%
    rename(barcode = variable) %>%
    filter(barcode %in% aliquot_barcode2cluster_df$barcode)
  
  ## map barcode to tumor subcluster
  cnv_state_mat.m$tumor_subcluster <- mapvalues(x = cnv_state_mat.m$barcode, from = aliquot_barcode2cluster_df$barcode, to = aliquot_barcode2cluster_df$ident)
  
  ## count number of cells with different cnv state per gene
  cnv_state_count <- cnv_state_mat.m %>%
    select(gene_symbol, tumor_subcluster, cna_state) %>%
    table() %>%
    as.data.frame() %>%
    filter(Freq > 0) %>%
    mutate(aliquot = aliquot_tmp)
  
  ## summarize the detected values for each gene
  cnv_nonna_count <- cnv_state_count %>%
    group_by(gene_symbol, tumor_subcluster) %>%
    summarize(num_cells_nonna = sum(Freq))
  
  ## merge counts
  cnv_state_count <- merge(cnv_state_count, cnv_nonna_count, by = c("gene_symbol", "tumor_subcluster"), all.x = T)
  
  ## combine with super table
  cnv_state_count_aliquots <- rbind(cnv_state_count, cnv_state_count_aliquots)
}
## calculate fraction
table2write <- cnv_state_count_aliquots
table2write$Fraction <- table2write$Freq/(table2write$num_cells_nonna)

## annotate each gene with expected cnv state and chromosome region
table2write <- merge(table2write, ccrcc_cna_genes_df, by = c("gene_symbol"))

## annotate cells with expected cnv state
table2write <- table2write %>%
  mutate(expected_cna = (gene_cna_type == "Loss" & cna_state %in% c(0.5, 0)) | (gene_cna_type == "Gain" & cna_state %in% c(1.5, 2)))

# write out CNA frequency table -------------------------------------------
write.table(x = table2write, file = paste0(dir_out, "Fraction_of_Malignant_Nephron_Epithelium_with_CNA_by_Gene.Per_Tumor_Subcluster.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)



