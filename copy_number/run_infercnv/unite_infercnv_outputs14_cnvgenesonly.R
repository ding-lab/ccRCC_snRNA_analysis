# Yige Wu @WashU Dec 2020
## plot InferCNV subcluster mode outputs onto UMAP for individual sample (all cells)

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

# input seurat object master list -----------------------------------------
## input CNV genes
ccrcc_cna_genes_df <- readxl::read_excel(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## input metadata
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")


# set run parameters ------------------------------------------------------
## specify CNV data
infercnv_run_id <- "20200305.v1"
dir_infercnv_all_runs <- "./Resources/snRNA_Processed_Data/InferCNV/outputs/"
dir_infercnv_output <- paste0(dir_infercnv_all_runs, "Individual.", infercnv_run_id, "/")
## set aliquot ids
aliquots2process <- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available]

# Loop: for each aliquot ---------
cnv_state_united_df <- NULL
# snRNA_aliquot_id_tmp <- "CPT0086820004"
for (snRNA_aliquot_id_tmp in aliquots2process) {
  if (snRNA_aliquot_id_tmp %in% c("CPT0000890002", "CPT0075170013")) {
    next()
  }
  ## get the readable aliquot id
  easy_id_tmp <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == snRNA_aliquot_id_tmp]
  
  ## input infercnv CNV state results
  obs_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  dim(obs_cnv_state_mat)
  dim(ref_cnv_state_mat)
  
  ## transform infercnv result wide data frame to long data frame
  cnv_state_df <- rbind(melt(obs_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))
  rm(obs_cnv_state_mat)
  rm(ref_cnv_state_mat)
  
  cnv_state_df <- cnv_state_df %>%
    filter(V1 %in% ccrcc_cna_genes_df$Gene_Symbol) %>%
    rename(gene_symbol = V1) %>%
    rename(barcode = variable) %>%
    rename(copy_state = value) %>%
    mutate(easy_id = easy_id_tmp) %>%
    mutate(aliquot = snRNA_aliquot_id_tmp) %>%
    select(easy_id, gene_symbol, barcode, copy_state, aliquot)
  
  ## unite
  cnv_state_united_df <- rbind(cnv_state_united_df, cnv_state_df)
}

# write outpu -------------------------------------------------------------
file2write <- paste0(dir_out, "infercnv.step14.outputs.", run_id, ".tsv")
write.table(file = file2write, x = cnv_state_united_df, row.names = F, sep = "\t", quote = F)
