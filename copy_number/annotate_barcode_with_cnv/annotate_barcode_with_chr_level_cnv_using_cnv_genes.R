# Yige Wu @WashU March 2020
## for annotating the barcode with cnv state per chromosome region using representative frequently altered genes (reported by TCGA)

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## set infercnv output directory
dir_infercnv_output <- "./Resources/snRNA_Processed_Data/InferCNV/outputs/"
# infercnv_run_id <- "Individual.20200305.v1"
infercnv_run_id <- "Individual.20200207.v1"
dir_infercnv_run <- paste0(dir_infercnv_output, infercnv_run_id, "/")
## get aliquots to process
aliquots2process <- list.files(dir_infercnv_run)
## chromosome regions from which the CNVs are annotated here
chr_regions2process <- unique(ccrcc_cna_genes_df$chr_region)
chr_regions2process <- as.vector(chr_regions2process)

# process by aliquot ------------------------------------------------------
cnv_state_per_cell_per_chr_region_aliquots <- NULL
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
  
  cnv_state_mat.m <- merge(cnv_state_mat.m, ccrcc_cna_genes_df, by = c("gene_symbol"), all.x = T)
  
  cnv_state_per_cell_per_chr_region_long <- cnv_state_mat.m %>%
    rename(cna_state_value = value) %>%
    rename(barcode = variable) %>%
    group_by(barcode, chr_region) %>%
    summarize(cnv_state = ifelse(all(cna_state_value == 1), "Neutral",
                                 ifelse(gene_cna_type == "Loss", ifelse(any(cna_state_value < 1), "Expected", "Other"),
                                        ifelse(any(cna_state_value > 1), "Expected", "Other"))))
  cnv_state_per_cell_per_chr_region_wide <- dcast(data = cnv_state_per_cell_per_chr_region_long, formula = barcode ~ chr_region, value.var = "cnv_state")
  chr_regions2add <- chr_regions2process[!(chr_regions2process %in% colnames(cnv_state_per_cell_per_chr_region_wide))]
  cnv_state_per_cell_per_chr_region_wide[,chr_regions2add] <- NA
  cnv_state_per_cell_per_chr_region_wide <- cnv_state_per_cell_per_chr_region_wide %>%
    mutate(aliquot = aliquot_tmp) %>%
    select(aliquot, barcode, chr_regions2process)
  cnv_state_per_cell_per_chr_region_wide %>% colnames()
  cnv_state_per_cell_per_chr_region_aliquots %>% colnames()
  cnv_state_per_cell_per_chr_region_aliquots <- rbind(cnv_state_per_cell_per_chr_region_wide, cnv_state_per_cell_per_chr_region_aliquots)
}


# write table -------------------------------------------------------------
write.table(x = cnv_state_per_cell_per_chr_region_aliquots,
            file = paste0(dir_out, infercnv_run_id, ".CNV_State_By_Chr_Region_By_Barcode.Using_Representative_Genes.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

