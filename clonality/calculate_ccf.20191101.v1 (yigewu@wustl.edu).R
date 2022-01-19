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
gene_cna_state_tab <- fread("Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/calculate_cna_freq/20191030.v1/TCGA_CNA_Genes_snCNA_Frequency.tsv", data.table = F)

# input CNA CCF by chromosome arm results ---------------------------
chr_amr_cna_state_tab <- fread("Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/calculate_cna_freq/20191030.v1/TCGA_CNA_Chr_Arm_snCNA_Frequency.tsv", data.table = F)

# set samples to process ----------------------------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# using VHL mutation or BAP1 or SETD2 or PBRM1 mutation to estimate tumor purity by WES sample -------------------------------------
mut_ccf_tab <- NULL
## create vectors to store the CCFs for different copy number states
tumor_perc_snRNA_aliquot_ids <- NULL
tumor_perc_genes <- NULL
tumor_perc_0xs <- NULL
tumor_perc_0.5xs <- NULL
tumor_perc_1xs <- NULL
tumor_perc_1.5xs <- NULL
tumor_perc_2xs <- NULL
tumor_perc_over_2xs <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  # snRNA_aliquot_id_tmp <- "CPT0001260013"
  
  ## choose which mutated gene will be used for tumpor purity estiamte
  ## choose the highest VAF among VHL, PBRM1, SETD2 and BAP1
  vaf_tmp <- vaf_tab %>%
    filter(Specimen.ID.snRNA == snRNA_aliquot_id_tmp) %>%
    select(SMGs[["CCRCC"]])
  
  gene0 <- colnames(vaf_tmp)[which.max(vaf_tmp)]
  gene0
  vaf0 <- max(vaf_tmp, na.rm = T)
  vaf0
  
  # get the cancer cell fraction of different copy number for the gene (for example VHL) --------
  tumor_perc_0x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 0 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
  tumor_perc_0.5x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 0.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
  tumor_perc_1.5x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 1.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
  tumor_perc_2x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 2 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
  tumor_perc_over_2x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 3 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
  if (length(tumor_perc_over_2x) == 0) {
    tumor_perc_over_2x <- 0
  } 
  tumor_perc_1x <- (1-tumor_perc_0x-tumor_perc_0.5x-tumor_perc_1.5x-tumor_perc_2x-tumor_perc_over_2x)
  
  a0 <- (1-tumor_perc_0x)
  
  if (gene0 != "KDM5C") {
    b0 <- (1*tumor_perc_0.5x + 2*tumor_perc_1x + 3*tumor_perc_1.5x + 4*tumor_perc_2x)
  } else {
    b0 <- (1*tumor_perc_1x + 1.5*tumor_perc_1.5x + 2*tumor_perc_2x)
  }
  
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
      tumor_perc_0.5x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 0.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_1.5x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 1.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_2x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 2 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_over_2x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 3 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      if (length(tumor_perc_over_2x) == 0) {
        tumor_perc_over_2x <- 0
      } 
      
      tumor_perc_1x <- (1-tumor_perc_0x-tumor_perc_0.5x-tumor_perc_1.5x-tumor_perc_2x-tumor_perc_over_2x)
      
      a <- (1-tumor_perc_0x)
      if (gene_tmp != "KDM5C") {
        b <- (1*tumor_perc_0.5x + 2*tumor_perc_1x + 3*tumor_perc_1.5x + 4*tumor_perc_2x)
      } else {
        b <- (1*tumor_perc_1x + 1.5*tumor_perc_1.5x + 2*tumor_perc_2x)
      }
      
      v <- as.numeric(vaf_tab %>%
                        filter(Specimen.ID.snRNA == snRNA_aliquot_id_tmp) %>%
                        select(gene_tmp))
      
      ccf <- (v*(b*p+2-2*p))/(p*a)
      ccf_mat[i,gene_tmp] <- ccf
      
      ## store the CCFs for different copy number states
      tumor_perc_snRNA_aliquot_ids <- c(tumor_perc_snRNA_aliquot_ids, snRNA_aliquot_id_tmp)
      tumor_perc_genes <- c(tumor_perc_genes, gene_tmp)
      tumor_perc_0xs <- c(tumor_perc_0xs, tumor_perc_0x)
      tumor_perc_0.5xs <- c(tumor_perc_0.5xs, tumor_perc_0.5x)
      tumor_perc_1xs <- c(tumor_perc_1xs, tumor_perc_1x)
      tumor_perc_1.5xs <- c(tumor_perc_1.5xs, tumor_perc_1.5x)
      tumor_perc_2xs <- c(tumor_perc_2xs, tumor_perc_2x)
      tumor_perc_over_2xs <- c(tumor_perc_over_2xs, tumor_perc_over_2x)
    }
  }
  mut_ccf_tab_tmp <- data.frame(snRNA_aliquot_id = snRNA_aliquot_id_tmp, gene0 = gene0, ccf0 = ccf0_to_test, tumor_purity = tumor_puritys, ccf_mat)
  mut_ccf_tab <- rbind(mut_ccf_tab_tmp, mut_ccf_tab)
}


# write mutation CCF table ------------------------------------------------
write.table(x = mut_ccf_tab, file = paste0(dir_out, "Somatic_Mutation_CCF_CNmut1.", run_id, ".tsv"), quote = F, row.names = F, col.names = T, sep = "\t")

# write CCFs for different copy number states for SMGs (easy to check what might be wrong after obtaining the mutation CCF) --------------------
tumor_perc_smg_cna_tab <- data.frame(snRNA_aliquot_id = tumor_perc_snRNA_aliquot_ids, gene_symbol = tumor_perc_genes, 
                                     tumor_perc_0x = tumor_perc_0xs, tumor_perc_0.5x = tumor_perc_0.5xs,
                                     tumor_perc_1x = tumor_perc_1xs,
                                     tumor_perc_1.5x = tumor_perc_1.5xs, tumor_perc_2x = tumor_perc_2xs,
                                     tumor_perc_over_2x = tumor_perc_over_2xs)
tumor_perc_smg_cna_tab <- unique(tumor_perc_smg_cna_tab)
write.table(x = tumor_perc_smg_cna_tab, file = paste0(dir_out, "SMG_CNA_CCF.", run_id, ".tsv"), quote = F, row.names = F, col.names = T, sep = "\t")


# pick the tumor purity estimated for WES closest to snRNA estimate tumor purity and corresponding CCFs ----------------
## for holding the final mutation CCF to present
mut_ccf_final_tab <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  snRNA_tumor_purity <- chr_amr_cna_state_tab %>%
    filter(snRNA_aliquot_id == snRNA_aliquot_id_tmp) %>%
    top_n(wt = perc_cna_in_all_cell, n = 1) %>%
    select(perc_cna_in_all_cell)
  
  snRNA_tumor_purity <- as.numeric(snRNA_tumor_purity)
  snRNA_tumor_purity <- min(snRNA_tumor_purity, 1)
  
  mut_ccf_final_tmp <- mut_ccf_tab %>%
    filter(snRNA_aliquot_id == snRNA_aliquot_id_tmp) %>%
    mutate(offset2snRNA_umor_purity = abs(tumor_purity - snRNA_tumor_purity)) %>%
    filter(tumor_purity <= 1) %>%
    filter(VHL <= 1) %>%
    filter(PBRM1 <= 1) %>%
    filter(SETD2 <= 1) %>%
    filter(KDM5C <= 1) %>%
    filter(PTEN <= 1) %>%
    filter(BAP1 <= 1) %>%
    filter(MTOR <= 1) %>%
    filter(TP53 <= 1) %>%
    top_n(wt = -offset2snRNA_umor_purity, n = 1)
  
  if (nrow(mut_ccf_final_tmp) == 0) {
    ## if there's no scenario that's possible, potentially the tumor_perc for copy number status is wrong
    ### pick the gene that have wrong CCF in most scenarios
    mut_ccf_mat2correct <- mut_ccf_tab %>%
      filter(snRNA_aliquot_id == snRNA_aliquot_id_tmp) %>%
      select(SMGs[["CCRCC"]])
    gene2correct <- colSums(mut_ccf_mat_tmp > 1)
    gene2correct <- names(gene2correct)[which.max(gene2correct)]
    gene2correct
    
    if (gene2correct %in% c("PBRM1")) {
      ## for this case I'll just  use the percerntage used for VHL
      
      vaf_tmp <- vaf_tab %>%
        filter(Specimen.ID.snRNA == snRNA_aliquot_id_tmp) %>%
        select(SMGs[["CCRCC"]])
      
      gene0 <- colnames(vaf_tmp)[which.max(vaf_tmp)]
      gene0
      vaf0 <- max(vaf_tmp, na.rm = T)
      vaf0
      
      # get the cancer cell fraction of different copy number for the gene (for example VHL) --------
      tumor_perc_0x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 0 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_0.5x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 0.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_1.5x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 1.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_2x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 2 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      tumor_perc_over_2x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene0 & gene_cna_state_tab$cna_state == 3 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
      if (length(tumor_perc_over_2x) == 0) {
        tumor_perc_over_2x <- 0
      } 
      tumor_perc_1x <- (1-tumor_perc_0x-tumor_perc_0.5x-tumor_perc_1.5x-tumor_perc_2x-tumor_perc_over_2x)
      
      a0 <- (1-tumor_perc_0x)
      
      if (gene0 != "KDM5C") {
        b0 <- (1*tumor_perc_0.5x + 2*tumor_perc_1x + 3*tumor_perc_1.5x + 4*tumor_perc_2x)
      } else {
        b0 <- (1*tumor_perc_1x + 1.5*tumor_perc_1.5x + 2*tumor_perc_2x)
      }
      
      
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
          if (gene_tmp != gene2correct) {
            tumor_perc_0x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 0 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
            tumor_perc_0.5x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 0.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
            tumor_perc_1.5x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 1.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
            tumor_perc_2x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 2 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
            tumor_perc_over_2x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == gene_tmp & gene_cna_state_tab$cna_state == 3 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
          } else {
            tumor_perc_0x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == "VHL" & gene_cna_state_tab$cna_state == 0 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
            tumor_perc_0.5x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == "VHL" & gene_cna_state_tab$cna_state == 0.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
            tumor_perc_1.5x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == "VHL" & gene_cna_state_tab$cna_state == 1.5 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
            tumor_perc_2x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == "VHL" & gene_cna_state_tab$cna_state == 2 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
            tumor_perc_over_2x <- as.numeric(gene_cna_state_tab$perc_cna_in_tumor_cell[gene_cna_state_tab$gene_symbol == "VHL" & gene_cna_state_tab$cna_state == 3 & gene_cna_state_tab$snRNA_aliquot_id == snRNA_aliquot_id_tmp])
          }
           if (length(tumor_perc_over_2x) == 0) {
            tumor_perc_over_2x <- 0
          } 
          
          tumor_perc_1x <- (1-tumor_perc_0x-tumor_perc_0.5x-tumor_perc_1.5x-tumor_perc_2x-tumor_perc_over_2x)
          
          a <- (1-tumor_perc_0x)
          if (gene_tmp != "KDM5C") {
            b <- (1*tumor_perc_0.5x + 2*tumor_perc_1x + 3*tumor_perc_1.5x + 4*tumor_perc_2x)
          } else {
            b <- (1*tumor_perc_1x + 1.5*tumor_perc_1.5x + 2*tumor_perc_2x)
          }
          
          v <- as.numeric(vaf_tab %>%
                            filter(Specimen.ID.snRNA == snRNA_aliquot_id_tmp) %>%
                            select(gene_tmp))
          
          ccf <- (v*(b*p+2-2*p))/(p*a)
          ccf_mat[i,gene_tmp] <- ccf
        }
      }
      mut_ccf_tab_tmp <- data.frame(snRNA_aliquot_id = snRNA_aliquot_id_tmp, gene0 = gene0, ccf0 = ccf0_to_test, tumor_purity = tumor_puritys, ccf_mat)
      
      mut_ccf_final_tmp <- mut_ccf_tab_tmp %>%
        mutate(offset2snRNA_umor_purity = abs(tumor_purity - snRNA_tumor_purity)) %>%
        filter(tumor_purity <= 1) %>%
        filter(VHL <= 1) %>%
        filter(PBRM1 <= 1) %>%
        filter(SETD2 <= 1) %>%
        filter(KDM5C <= 1.00) %>%
        filter(PTEN <= 1) %>%
        filter(BAP1 <= 1) %>%
        filter(MTOR <= 1) %>%
        filter(TP53 <= 1) %>%
        top_n(wt = -offset2snRNA_umor_purity, n = 1)
      
    } else {
      stop("I don't know what to do!")
    }
    
  }
    
  mut_ccf_final_tab <- rbind(mut_ccf_final_tmp, mut_ccf_final_tab)
}
write.table(x = mut_ccf_final_tab, file = paste0(dir_out, "Somatic_Mutation_CCF_CNmut1.Closest2snRNA_TumorPurity.", run_id, ".tsv"), quote = F, row.names = F, col.names = T, sep = "\t")


# Trasnform the Chr Arm CNA CCF to Matrix Like ----------------------------
chr_amr_cna_state_tab <- chr_amr_cna_state_tab %>%
  mutate(cna_text = paste0(chr_region, "_", gene_cna_type)) %>%
  mutate(snRNA_tumor_purity = number_tumor_cells/num_all_cell)
chr_amr_cna_state_mat <- dcast(data = chr_amr_cna_state_tab, formula = snRNA_aliquot_id ~ cna_text, value.var = "perc_cna_in_tumor_cell")
chr_amr_cna_state_mat[,-1][chr_amr_cna_state_mat[,-1] > 1] <- 1
snRNA_tumor_purity <- chr_amr_cna_state_tab %>%
  group_by(snRNA_aliquot_id) %>%
  top_n(wt = perc_cna_in_all_cell, n = 1) %>%
  select(snRNA_aliquot_id, perc_cna_in_all_cell) %>%
  rename(snRNA_tumor_purity = perc_cna_in_all_cell)

chr_amr_cna_state_mat <- chr_amr_cna_state_mat %>%
  mutate(snRNA_tumor_purity = as.numeric(mapvalues(x = chr_amr_cna_state_mat$snRNA_aliquot_id, from = snRNA_tumor_purity$snRNA_aliquot_id, to = snRNA_tumor_purity$snRNA_tumor_purity)))
write.table(x = chr_amr_cna_state_mat, file = paste0(dir_out, "TCGA_CNA_Chr_Arm_snRNA_CCF_Matrix.", run_id, ".tsv"), quote = F, row.names = F, col.names = T, sep = "\t")

# Transform the Mutation CCF to be ready for merging---------------------------------------------------------------
smg_to_keep <- SMGs[["CCRCC"]]
smg_to_keep <- smg_to_keep[colSums(mut_ccf_final_tab[,SMGs[["CCRCC"]]]) > 0]
mut_ccf_tab2merge <- mut_ccf_final_tab %>%
  select("snRNA_aliquot_id", smg_to_keep, "tumor_purity") %>%
  rename(bulk_WES_tumor_purity = tumor_purity)

# merge Mutation CCF with CNA CCF -----------------------------------------
ccf_tab <- data.frame(Aliquot = chr_amr_cna_state_mat$snRNA_aliquot_id)
## add Case Id
ccf_tab$Case <- as.vector(mapvalues(x = ccf_tab$Aliquot, from = meta_tab$Specimen.ID.snRNA, to = meta_tab$Case.ID))
ccf_tab$Case[ccf_tab$Aliquot == "CPT0025890002"] <- "C3N-00733"
ccf_tab$Case[ccf_tab$Aliquot == "CPT0075140002"] <- "C3N-01200"

## input clinical table
clinial_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/generate_clinical_table/20191017.v1/snRNA_ccRCC_Clinicl_Table.20191017.v1.tsv", data.table = F)

## add Gender
ccf_tab$Gender <- mapvalues(x = ccf_tab$Case, from = clinial_tab$Case, to = clinial_tab$Gender)
ccf_tab <- ccf_tab %>%
  select(Case, Gender, Aliquot)
ccf_tab <- merge(ccf_tab, chr_amr_cna_state_mat, by.x = "Aliquot", by.y = c("snRNA_aliquot_id"), all.x = T)
ccf_tab <- merge(ccf_tab, mut_ccf_tab2merge, by.x = "Aliquot", by.y = c("snRNA_aliquot_id"), all.x = T)
# write out results -------------------------------------------------------
write.table(x = ccf_tab, file = paste0(dir_out, "Somatic_Mutation_Chr_Arm_CNA_CCF_CNmut1.Closest2snRNA_TumorPurity.", run_id, ".tsv"), quote = F, row.names = F, col.names = T, sep = "\t")


