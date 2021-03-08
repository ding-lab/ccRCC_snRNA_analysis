# Yige Wu @ WashU 2019 Mar
## extract VAF info for MAF files for prospective BRCA and CO

# extract VAF info from vcf files -----------------------------------------
## MGI: gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/somatic/BRCA/
## MGI: gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/somatic/CO/
## downloaded to Box: https://wustl.app.box.com/folder/69912658362

# source ------------------------------------------------------------------
setwd(dir = "~/Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")


# download vcfs from MGI --------------------------------------------------
## cat /Users/yigewu/Box Sync/Ding_Lab/Projects_Current/ccRCC_drug/resources/Somatic_Mutation/scRNA_RCC_test/worklog

# set variables -----------------------------------------------------------
gene <- "TP53"


# input the patient id conversion table -----------------------------------
partID_map_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/SampleMatch_For_ExomeData.121517.tsv", data.table = F)

maf_gene_new <- NULL
maf <- loadMaf(cancer = cancer_tmp, maf_files = maf_files)

maf_gene <- maf %>%
  filter(Hugo_Symbol == gene) %>%
  mutate(variant_id = paste0(Chromosome, "_", Start_Position, "_", End_Position)) %>%
  mutate(partID = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", n = 2)[,1])
row_tmp <- which(maf_gene$Variant_Type == "DEL")
maf_gene$variant_id[row_tmp] <- paste0(maf_gene$Chromosome[row_tmp], "_", as.numeric(maf_gene$Start_Position[row_tmp]) - 1, "_", maf_gene$End_Position[row_tmp])

for (partID_tmp in "RCC_A1") {
  maf_gene_partID <- maf_gene %>%
    filter(partID == partID_tmp)
  for (variant_id_tmp in unique(maf_gene_partID$variant_id)) {
    maf_gene_partID_variant <- maf_gene_partID %>%
      filter(variant_id == variant_id_tmp)
    
    vcf_path_tmp <- paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/per_sample/gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/somatic/", 
                           cancer_tmp, "/", partID_tmp, "/merged.vcf")
    if (!file.exists(vcf_path_tmp)) {
      partID_tmp_new <- partID_map_tab$SampleID_In_ExomeData[partID_map_tab$Sample_Name == partID_tmp]
      partID_tmp_old <- partID_tmp
      partID_tmp <- partID_tmp_new
      vcf_path_tmp <- paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/per_sample/gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/somatic/", 
                             cancer_tmp, "/", partID_tmp, "/merged.vcf")
    }
    
    vcf_file_tmp <- fread(input = vcf_path_tmp, data.table = F)
    vcf_file_tmp <- fread(input = "./Ding_Lab/Projects_Current/ccRCC_drug/resources/Somatic_Mutation/scRNA_RCC_test/RCC_A1.vep_filtered.vcf", data.table = F, fill = T, )
    
    vcf_file_tmp <- vcf_file_tmp %>%
      mutate(End_Position = POS) 
    
    row_tmp <- (nchar(vcf_file_tmp$REF) > 1)
    vcf_file_tmp$End_Position[row_tmp] <- (as.numeric(vcf_file_tmp$POS[row_tmp]) + nchar(vcf_file_tmp$REF[row_tmp]) - 1)
    vcf_file_tmp$Chromosome <- vcf_file_tmp$`#CHROM`
    vcf_file_tmp <- vcf_file_tmp %>%
      mutate(variant_id = paste0(Chromosome, "_", POS, "_", End_Position))
    
    vcf_gene_tmp <- vcf_file_tmp %>%
      filter(variant_id == variant_id_tmp)
    vaf_tmp <- NA
    if (nrow(vcf_gene_tmp) > 0) {
      if (grepl(pattern = "AU", vcf_gene_tmp$FORMAT)) {
        ## streka-varscan format
        t_alt_tmp <- vcf_gene_tmp$ALT
        t_format_tmp <- str_split(string = vcf_gene_tmp$FORMAT, pattern = ":")[[1]]
        t_count_tmp <- str_split(string = vcf_gene_tmp$TUMOR, pattern = ":")[[1]]
        t_alt_count_tmp <- t_count_tmp[which(t_format_tmp == paste0(t_alt_tmp, "U"))]
        t_alt_count_tmp <- as.numeric(str_split(string = t_alt_count_tmp, pattern = ",")[[1]][1])
        t_dp_tmp <- as.numeric(t_count_tmp[which(t_format_tmp == "DP")])
        vaf_tmp <- (t_alt_count_tmp/t_dp_tmp)
      } else if (grepl(pattern = "AD", vcf_gene_tmp$FORMAT)) {
        t_format_tmp <- str_split(string = vcf_gene_tmp$FORMAT, pattern = ":")[[1]]
        t_count_tmp <- str_split(string = vcf_gene_tmp$TUMOR, pattern = ":")[[1]]
        t_alt_count_tmp <- as.numeric(t_count_tmp[which(t_format_tmp == "AD")])
        t_dp_tmp <- as.numeric(t_count_tmp[which(t_format_tmp == "DP")])
        vaf_tmp <- (t_alt_count_tmp/t_dp_tmp)
      } 
    }
    
    maf_gene_partID_variant$VAF <- vaf_tmp
    maf_gene_new <- rbind(maf_gene_new, maf_gene_partID_variant)
    
  }
}