# map_infercnv_state2category <- function(copy_state) {
#   cnv_cat <- vector(mode = "character", length = length(copy_state))
#   cnv_cat[is.na(copy_state)] <- "Not Available"
#   cnv_cat[copy_state == 1] <- "Neutral"
#   cnv_cat[copy_state == 0] <- "Complete Loss"
#   cnv_cat[copy_state == 0.5] <- "Loss of one copy"
#   cnv_cat[copy_state == 1.5] <- "Addition of one copy"
#   cnv_cat[copy_state == 2] <- "Addition of two copies"
#   cnv_cat[copy_state == 3] <- "Addition > two copies"
#   return(cnv_cat)
# }
# 
# 
# copy_number_colors <-  c("Complete Loss" = PuBu_colors[9],
#                    "Loss of one copy" = PuBu_colors[5],
#                    "Neutral" = PuBu_colors[3],
#                    "Addition of one copy" = PuRd_colors[5], 
#                    "Addition of two copies" = PuRd_colors[7],
#                    "Addition > two copies" = PuRd_colors[9],
#                    "Not Available" = "grey50")


# load omics data ---------------------------------------------------------
loadMaf <- function() {
  maf <- fread(input = paste0(dir_base, "Resources/Bulk_Processed_Data/Somatic_Variants/ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf"), data.table = F, fill=TRUE)
  maf <- data.frame(maf)
  print(paste0("MAF has ", nrow(maf), " lines\n"))
  return(maf)
}

loadCNAstatus <- function() {
  ## input CNA values
  cna <- fread(input = paste0(dir_base, "Resources/Bulk_Processed_Data/WGS_CNV_Somatic/Michigan/somatic_CNA.CCRCC.partID.txt"), data.table = F)
  ## reformat
  cna_head <- cna$gene
  cna_mat <- cna[, colnames(cna)[!(colnames(cna) %in% "gene")]]
  ## create new matrix for the CNV status
  cna_status <- matrix(data = "neutral", nrow = nrow(cna_mat), ncol = ncol(cna_mat))
  cna_status[cna_mat > 0.1] <- "amplification"
  cna_status[cna_mat < -0.1] <- "deletion"
  cna_status <- data.frame(cbind(cna$gene, cna_status))
  colnames(cna_status) <- colnames(cna)
  return(cna_status)
}
# 
# loadRNA <- function() {
#   rna_fn <- "ccRcc_RNA_rpkm_Mich_formatted_tumor.csv"
#   rna_tab <- fread(input = paste0(dir2dinglab_projects, "TP53_shared_data/resources/rna/", rna_fn), data.table = F)
#   
#   colnames(rna_tab)[1] <- "gene"
#   rna_mat <- as.matrix(rna_tab[,-1])
#   rna_mat_log2 <- log2(rna_mat+1)
#   rna_tab <- data.frame(gene = rna_tab$gene)
#   rna_tab <- cbind(rna_tab, as.data.frame(rna_mat_log2))
#   return(rna_tab) 
# }
# 
loadParseProteomicsData <- function(expression_type, sample_type) {
  ## expresson_type: PRO or PHO (phosphosite level) or PHO_collapsed (protein level)
  if (expression_type == "protein") {
    path_protein <- paste0("./Resources/Bulk_Processed_Data/Protein/CCRCC_PRO_tumor_PGDAC_MD_MAD_partID.txt")
    exp_data <- fread(input = path_protein, data.table = F)
  }
  return(exp_data)
}
