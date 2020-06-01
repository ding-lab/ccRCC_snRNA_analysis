# Yige Wu @WashU May 2020

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
## input DEGs
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findallmarkers_wilcox_tumorcells_by_manualsubcluster/20200427.v1/Tumormanualsubcluster.FindAllMarkers.Wilcox.Minpct0.1.Logfc0.25.tsv")
nrow(markers_df)
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## input known cnv genes
known_cnv_genes_df <- readxl::read_xlsx(path = "./Resources/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## input ppi table
ppi_pair_df <- fread(data.table = F, input = "./Resources/Databases/Protein_Protein_Interactions/protein_pair_table_v2.txt")
## load CNV fraction in tumor cells
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20200512.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20200512.v1.tsv", data.table = F)
table(cnv_3state_count_aliquots$tumor_subcluster)

# get genes that are altered for each sampe --------------------------------------------
## add aliquot.wu
cnv_3state_count_aliquots$aliquot.wu <- mapvalues(x = cnv_3state_count_aliquots$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
cna_frac_df <- cnv_3state_count_aliquots
cna_frac_df$case <- mapvalues(x = cna_frac_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))

## add cytoband and expected cna type
cna_frac_df$gene_cytoband <- mapvalues(x = cna_frac_df$gene_symbol, from = known_cnv_genes_df$Gene_Symbol, to = as.vector(known_cnv_genes_df$Cytoband))
cna_frac_df$gene_expected_state <- mapvalues(x = cna_frac_df$gene_symbol, from = known_cnv_genes_df$Gene_Symbol, to = as.vector(known_cnv_genes_df$CNV_Type))
cna_frac_df <- cna_frac_df %>%
  mutate(id_aliquot_cluster = paste0(aliquot.wu, "_C", (tumor_subcluster + 1)))

## get the data with only expected CNV state
cna_frac_filtered_df <- cna_frac_df %>%
  filter(gene_expected_state == cna_3state) 

## filter genes
genes_filtered <- unique(cna_frac_filtered_df$gene_symbol[cna_frac_filtered_df$Fraction > 0.1])
cna_frac_filtered_df <- cna_frac_filtered_df %>%
  filter(gene_symbol %in% genes_filtered) %>%
  mutate(Fraction_Range = ifelse(Fraction < 0.5, "<=50%", ">50%")) %>%
  mutate(aliquot.wu = str_split_fixed(string = id_aliquot_cluster, pattern = "_", n = 2)[,1]) %>%
  mutate(name_tumorsubcluster = str_split_fixed(string = id_aliquot_cluster, pattern = "_", n = 2)[,2]) %>%
  mutate(Data_detected = T) %>%
  dplyr::select(id_aliquot_cluster, gene_symbol, Fraction, name_tumorsubcluster, cna_3state, Fraction_Range, Data_detected, aliquot.wu, gene_cytoband)

# format marker table -----------------------------------------------------
## filter markers
markers_df <- markers_df %>%
  filter(p_val_adj < 0.05)
## add cluster name 
markers_df$id_aliquot_wu <- mapvalues(x = markers_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
markers_df <- markers_df %>%
  mutate(name_cluster = paste0(id_aliquot_wu, "_C", (cluster + 1)))

# format the proteins known to interact with proteins from cna genes -----------------------------------
ppi_pair_cnvgenes_df <- merge(ppi_pair_df, known_cnv_genes_df, by.x = c("GENE"), by.y = c("Gene_Symbol"), all.y = T)

# filter degs by cna interactome ------------------------------------------
# markers_filtered_df <- markers_df %>%
#   filter(gene %in% ppi_pair_cnvgenes_df$SUB_GENE)
markers_filtered_df <- merge(markers_df, ppi_pair_cnvgenes_df %>%
                               dplyr::rename(cna_gene_symbol = GENE), by.x = c("gene"), by.y = c("SUB_GENE"))
## filter by genes altered in each sample
markers_filtered_df <- merge(markers_filtered_df, 
                             cna_frac_filtered_df %>%
                               dplyr::select(id_aliquot_cluster, gene_symbol), by.x = c("name_cluster", "cna_gene_symbol"), by.y = c("id_aliquot_cluster", "gene_symbol"))
## divide deg into up and down-regulated
markers_filtered_df <- markers_filtered_df %>%
  mutate(deg_direction = ifelse(avg_logFC > 0, "up", "down"))
count_deg_by_ppi_df <- markers_filtered_df %>%
  dplyr::select(name_cluster, deg_direction, cna_gene_symbol) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0)

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "tumormanualsubcluster.", "deg_count.", "cnv_genes_interactome.", run_id, ".tsv")
write.table(x = count_deg_by_ppi_df, file = file2write, quote = F, sep = "\t", row.names = F)

