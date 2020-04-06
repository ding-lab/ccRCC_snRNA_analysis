# Yige Wu @WashU Oct 2019
## merge CNA frequency data (inferCNV) with mutation frequency data (WES VAF)

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input the InferCNV processed CNA frequency table ------------------------
infercnv_freq_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/heatmap_infercnv/20191028.v1/TCGA_CNA_Genes_snCNA_Frequency.tsv", data.table = F)

# estimate 3p, 5q, 14q CNA frequency by pathogenic target genes -----------
infercnv_freq_tab.collapse.sn_cna_cat <- infercnv_freq_tab %>%
  group_by(aliquot, cna_text, chr_region) %>%
  summarize(perc_cna_in_all_cell.collapsed = sum(perc_cna_in_all_cell, na.rm = T))

infercnv_freq_tab.chr_region <- infercnv_freq_tab.collapse.sn_cna_cat %>%
  group_by(aliquot, chr_region) %>%
  summarize(perc_cna_in_all_cell.top.chr_region = max(perc_cna_in_all_cell.collapsed))
  
# transform CNA frequency table to matrix like ----------------------------
infercnv_freq_tab.chr_region.mat <- dcast(data = infercnv_freq_tab.chr_region, formula = aliquot ~ chr_region, value.var = c("perc_cna_in_all_cell.top.chr_region"))
write.table(x = infercnv_freq_tab.chr_region.mat, file = paste0(dir_out, "TCGA_CNA_Chr_Regions_snCNA_Frequency.tsv"), quote = F, sep = "\t", row.names = F)

# input mutation VAF table ------------------------------------------------
vaf_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/mutation/generate_bulk_mutation_table/20191024.v1/snRNA_ccRCC_Mutation_VAF_Table.20191024.v1.csv", data.table = F)


# change VAF to frequency -------------------------------------------------
vaf_inferred_freq_tab <- vaf_tab[, SMGs[["CCRCC"]]]
vaf_inferred_freq_tab <- 2*vaf_inferred_freq_tab
vaf_inferred_freq_tab <- cbind(vaf_tab[, c("Case", "Aliquot")], vaf_inferred_freq_tab)

# load meta data ----------------------------------------------------------
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)

# merge CNA with VAF ------------------------------------------------------
freq_sup_tab <- merge(infercnv_freq_tab.chr_region.mat, meta_tab, by.x = c("aliquot"), by.y = c("Specimen.ID.snRNA"), all.x = T)
freq_sup_tab$Case.ID[freq_sup_tab$aliquot == "CPT0075140002"] <- "C3N-01200"
freq_sup_tab$Case.ID[freq_sup_tab$aliquot == "CPT0025890002"] <- "C3N-00733"

freq_sup_tab <- merge(freq_sup_tab, vaf_tab, by.x = c("Case.ID"), by.y = c("Case"), all.x = T)
write.table(x = freq_sup_tab, file = paste0(dir_out, "TCGA_CNA_Chr_Regions_snCNA_Frequency_with_VAF.tsv"), quote = F, sep = "\t", row.names = F)

# make simplified CNA frequency with VAFtable ---------------------------------------------------
sim_freq_tab <- freq_sup_tab %>%
  select(c("Case.ID", "aliquot", "3p", "5q", "14q", SMGs[["CCRCC"]]))
write.table(x = sim_freq_tab, file = paste0(dir_out, "Simplied_TCGA_CNA_Chr_Regions_snCNA_Frequency_with_VAF.tsv"), quote = F, sep = "\t", row.names = F)

# merge CNA with VAF*2 ------------------------------------------------------
freq_sup_tab <- merge(infercnv_freq_tab.chr_region.mat, meta_tab, by.x = c("aliquot"), by.y = c("Specimen.ID.snRNA"), all.x = T)
freq_sup_tab$Case.ID[freq_sup_tab$aliquot == "CPT0075140002"] <- "C3N-01200"
freq_sup_tab$Case.ID[freq_sup_tab$aliquot == "CPT0025890002"] <- "C3N-00733"

freq_sup_tab <- merge(freq_sup_tab, vaf_inferred_freq_tab, by.x = c("Case.ID"), by.y = c("Case"), all.x = T)
write.table(x = freq_sup_tab, file = paste0(dir_out, "TCGA_CNA_Chr_Regions_snCNA_Frequency_with_VAF_inferred_Mut_Freq.tsv"), quote = F, sep = "\t", row.names = F)


# make simplified frequency table ---------------------------------------------------
sim_freq_tab <- freq_sup_tab %>%
  select(c("Case.ID", "aliquot", "3p", "5q", "14q", SMGs[["CCRCC"]]))
write.table(x = sim_freq_tab, file = paste0(dir_out, "Simplied_TCGA_CNA_Chr_Regions_snCNA_Frequency_with_VAF_inferred_Mut_Freq.tsv"), quote = F, sep = "\t", row.names = F)



