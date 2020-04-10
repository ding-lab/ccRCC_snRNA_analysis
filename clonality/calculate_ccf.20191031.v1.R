# Yige Wu @WashU Oct 2019
## for calculating cancer cell fraction of CNAs and Mutations


# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input VAF table ---------------------------------------------------------
vaf_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/mutation/generate_bulk_mutation_table/20191024.v1/snRNA_ccRCC_Mutation_VAF_Table.20191024.v1.csv", data.table = F)

# load meta data ----------------------------------------------------------
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)

# add snRNA_aliquot_id to vaf table ---------------------------------------
vaf_tab <- merge(vaf_tab, meta_tab, by.x = c("Aliquot"), by.y = c("Specimen.ID.bulk"), all.x = T)

# input CNA frequency by gene ---------------------------------------------
# gene_cna_state_tab <- 

# set samples to process ----------------------------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# using VHL mutation or BAP1 or SETD2 or PBRM1 mutation to estimate tumor purity by WES sample -------------------------------------
ccf_tab <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  # snRNA_aliquot_id_tmp <- "CPT0001260013"
  
  ## choose which mutated gene will be used for tumpor purity estiamte
  ## choose the highest VAF among VHL, PBRM1, SETD2 and BAP1
  vaf_tmp <- vaf_tab %>%
    filter(Specimen.ID.snRNA == snRNA_aliquot_id_tmp) %>%
    select(VHL, PBRM1, BAP1, SETD2)
  gene0 <- colnames(vaf_tmp)[which.max(vaf_tmp)]
  gene0
  vaf0 <- max(vaf_tmp, na.rm = T)
  vaf0
  
  # get the cancer cell fraction of different copy number for the gene (for example VHL) --------
  tumor_perc_0x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 0 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
  tumor_perc_1x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 0.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
  tumor_perc_3x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 1.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
  tumor_perc_4x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 2 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
  tumor_perc_over_4x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 3 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
  if (length(tumor_perc_over_4x) != 0) {
    tumor_perc_2x <- (1-tumor_perc_0x-tumor_perc_1x-tumor_perc_3x-tumor_perc_4x-tumor_perc_over_4x)
  } else {
    tumor_perc_2x <- (1-tumor_perc_0x-tumor_perc_1x-tumor_perc_3x-tumor_perc_4x)
  }
  
  a0 <- (1-tumor_perc_0x)
  
  b0 <- (1*tumor_perc_1x + 2*tumor_perc_2x + 3*tumor_perc_3x + 4*tumor_perc_4x)
  
  ccf0_to_test <- seq(from = 1, to = 0, by = -0.05)
  ## create a vector to hold the estimated tumor purity
  tumor_puritys <- NULL
  ## create a matrix to hold the CCF for the rest of the mutated gene
  ccf_mat <- matrix(data = 0, nrow = length(ccf0_to_test), ncol = length(SMGs[["CCRCC"]]))
  colnames(ccf_mat) <- SMGs[["CCRCC"]]
  for (i in 1:length(ccf0_to_test)) {
    ccf0 <- ccf0_to_test[i]
    
    # calculate purity according to the assumed ccf for the gene (for exapmle VHL--------------------------------------------------------
    p <- 2/((a0*ccf0)/vaf0 + 2 - b0)
    tumor_puritys <- c(tumor_puritys, p) 
    
    genes2test <- vaf_tab %>%
      filter(Specimen.ID.snRNA == snRNA_aliquot_id_tmp)
    genes2test <- colnames(genes2test[1, !is.na(genes2test)])
    genes2test <- intersect(genes2test, SMGs[["CCRCC"]])
    
    for (gene_tmp in genes2test) {

      tumor_perc_0x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 0 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_1x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 0.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_3x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 1.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_4x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 2 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_over_4x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 3 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      if (length(tumor_perc_over_4x) != 0) {
        tumor_perc_2x <- (1-tumor_perc_0x-tumor_perc_1x-tumor_perc_3x-tumor_perc_4x-tumor_perc_over_4x)
      } else {
        tumor_perc_2x <- (1-tumor_perc_0x-tumor_perc_1x-tumor_perc_3x-tumor_perc_4x)
      }

      a <- (1-tumor_perc_0x)
      if (gene_tmp != "KDM5C") {
        b <- (1*tumor_perc_1x + 2*tumor_perc_2x + 3*tumor_perc_3x + 4*tumor_perc_4x)
      } else {
        b <- (1*tumor_perc_2x + 1.5*tumor_perc_3x + 2*tumor_perc_4x)
      }
      
      v <- as.numeric(vaf_tab %>%
                        filter(Specimen.ID.snRNA == snRNA_aliquot_id_tmp) %>%
                        select(gene_tmp))
      
      ccf <- (v*(b*p+2-2*p))/(p*a)
      ccf_mat[i,gene_tmp] <- ccf
    }
  }
  ccf_tab_tmp <- data.frame(snRNA_aliquot_id = snRNA_aliquot_id_tmp, gene0 = gene2use, ccf0 = ccf0_to_test, tumor_purity = tumor_puritys, ccf_mat)
  ccf_tab <- rbind(ccf_tab_tmp, ccf_tab)
}
write.table(x = ccf_tab, file = paste0(dir_out, "CCF_CNmut1.", run_id, ".tsv"), quote = F, row.names = F, col.names = T, sep = "\t")

# estimate CCF for somatic mutations using estimated tumor purity --------------------------------------

# input CNA CCF by chromosome arm results ---------------------------


# merge Mutation CCF with CNA CCF -----------------------------------------


# write out results -------------------------------------------------------


