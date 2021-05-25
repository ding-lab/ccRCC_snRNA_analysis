# Yige Wu @WashU May 2021

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
dap_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Differential_Peaks/20210510/da_up_peaks_Tumor_vs_PT.annotated.20210510.tsv")
## input the da peaks with hg19 coordinates
dap_hg19_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/enhancer/convert_up_dap_to_hg19/20210511.v1/Up_DAP.NonPromoter.withHg19.tsv")
## input enhancer prediction
enhancer_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Differential_Peaks/20210511/243670965_enh_target.txt", col.names = c("chr", "start.hg19", "end.hg19", "Target_Genes", "Cell_resource"))
enhancer_df2 <- readxl::read_excel(path = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Differential_Peaks/20210511/243670965_enh_target.xlsx")

# preprocess the enhancer prediction --------------------------------------
enhancer_df <- enhancer_df %>%
  mutate(peak.hg19 = paste0(chr, ":", start.hg19, "-", end.hg19))
enhancer2genes_df <- merge(x = enhancer_df, 
                     y = enhancer_df2 %>%
                       dplyr::rename(DAP_Type = Type) %>%
                       dplyr::select(Peaks, DAP_Type), by.x = c("peak.hg19"), by.y = c("Peaks"), all.x = T)
idx_rep <- sapply(X = enhancer2genes_df$Target_Genes, FUN = function(x) {
  genes_vec <- str_split(string = x, pattern = ";")[[1]]
  len_genes <- length(genes_vec)
  return(len_genes)
})
genes_uniq <- sapply(X = enhancer2genes_df$Target_Genes, FUN = function(x) {
  genes_vec <- str_split(string = x, pattern = ";")[[1]]
  return(genes_vec)
})
enhancer2gene_df <- enhancer2genes_df[rep(1:nrow(enhancer2genes_df), idx_rep),]
enhancer2gene_df$Target_Gene <- unlist(genes_uniq)

# merge -------------------------------------------------------------------
dap2gene_df <- merge(x = dap_df, y = dap_hg19_df, by = intersect(colnames(dap_df), colnames(dap_hg19_df)), all.x = T)
dap2gene_df <- merge(x = dap2gene_df %>%
                         mutate(peak.hg19 = paste0(chr, ":", start.hg19, "-", end.hg19)),
                       y = enhancer2gene_df, by = c("peak.hg19", "chr", "start.hg19", "end.hg19"), all.x = T)
dap2gene_sim_df <- dap2gene_df %>%
  filter(Count_up >= 5) %>%
  mutate(Peak_Type = ifelse(Type == "Promoter", "Promoter",
                                  ifelse(!is.na(DAP_Type) & DAP_Type == "Enhancer", "Enhancer", Type))) %>%
  mutate(Peak_Type = ifelse(grepl(x =  Peak_Type, pattern = "Intron"), "Intron",
                            ifelse(grepl(x = Peak_Type, pattern = "Exon"), "Exon",
                                   ifelse(grepl(x = Peak_Type, pattern = "Downstream"), "Downstream", Peak_Type)))) %>%
  mutate(Gene = ifelse(Peak_Type == "Enhancer", Target_Gene, Gene)) %>%
  rename(Peak = peak) %>%
  select(Peak, Gene, Peak_Type) %>%
  unique()
table(dap2gene_sim_df$Peak_Type)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DAPeak2Gene.",  run_id, ".tsv")
write.table(file = file2write, x = dap2gene_df, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "DAPeak2Gene.Count5.",  run_id, ".tsv")
write.table(file = file2write, x = dap2gene_sim_df, quote = F, sep = "\t", row.names = F)
