# Yige Wu @WashU March 2020
## for annotating the barcode with cnv state per chromosome region using representative frequently altered genes (reported by TCGA)

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
## get aliquots to process
files2process <- list.files(dir_infercnv_output, recursive = T)
files2process <- files2process[grepl(pattern = "CPT", x = files2process) & grepl(pattern = "Individual.20200305.v1|run.20210805", x = files2process)]
files2process <- files2process[grepl(pattern = "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt", x = files2process)]
files2process
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
genes2process <- knowncnvgenes_df$Gene_Symbol

# process by aliquot ------------------------------------------------------
cnv_state_bycell_bygene_df <- NULL
for (file_tmp in files2process) {
# for (aliquot_tmp in aliquots2process) {
  ## input infercnv cnv state outputs
  cnv_state_obs_mat <- fread(input = paste0(dir_infercnv_output, file_tmp), data.table = F)
  cnv_state_ref_mat <- fread(input = paste0(dir_infercnv_output, gsub(x = file_tmp, pattern = "observations", replacement = "references")), data.table = F)
  aliquot_tmp <- str_split(string = file_tmp, pattern = "\\/")[[1]][2]
  ## filter cnv results to only selected genes
  cnv_state_obs_mat <- cnv_state_obs_mat %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% genes2process)
  cnv_state_ref_mat <- cnv_state_ref_mat %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% genes2process)
  
  ## melt wide data frame to long data frame
  cnv_state_obs_mat.m <- melt(cnv_state_obs_mat, id.vars = c("gene_symbol"))
  cnv_state_ref_mat.m <- melt(cnv_state_ref_mat, id.vars = c("gene_symbol"))
  
  ## combine the cnv results
  cnv_state_mat.m <- rbind(cnv_state_obs_mat.m, cnv_state_ref_mat.m)
  cnv_state_mat.m <- cnv_state_mat.m %>%
    rename(barcode_individual = variable) %>%
    rename(cna_value = value) %>%
    mutate(gene_expected_cna_state = mapvalues(x = gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = knowncnvgenes_df$CNV_Type)) %>%
    mutate(cna_state = ifelse(cna_value == 1, "Neutral",
                              ifelse(cna_value < 1, "Loss", "Gain"))) %>%
    mutate(id_aliquot = aliquot_tmp)
  cnv_state_bycell_bygene_df <- rbind(cnv_state_bycell_bygene_df, cnv_state_mat.m)
}

# write table -------------------------------------------------------------
write.table(x = cnv_state_bycell_bygene_df,
            file = paste0(dir_out, "CNV_State_By_Gene_By_Barcode.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

