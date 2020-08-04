# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input pseudotime degs
deg_by_pseudime_df <- fread(data.table = F, input = "./Resources/Analysis_Results/monocle/deg/unite_pseudotime_degs_across_3samples/20200728.v1/PT_and_TumorCells_ByCellType.Pseudotim_DE_genes.Combine3Samples.20200728.v1.tsv")
## input degs between tumor and NATs in bulk data
deg_bulk_rna_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/compare_bulk_mrna_tumor_vs_normal/20200728.v1/Bulk_mRNA_Tumor_vs_Normal.Wilcox.20200728.v1.tsv")
deg_bulk_pro_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/compare_bulk_protein_tumor_vs_normal/20200728.v1/Bulk_Protein_Tumor_vs_Normal.Wilcox.20200728.v1.tsv")
## input seurat cell type degs (tumor vs normal PT)
deg_by_celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/examine_degs_roc_for_tumor_vs_normalpt/20200728.v1/findmarkers_roc.tumor_vs_normalpt.20200728.v1.tsv")
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv", data.table = F)

# filter pseudotime degs by use_for_ordering, qval and the presence across samples -----------------------------------
deg_by_pseudime_filtered_df <- deg_by_pseudime_df %>%
  filter(use_for_ordering) %>%
  filter(qval<0.1)
deg_count_df <- deg_by_pseudime_filtered_df %>%
  select(gene_short_name) %>%
  table() %>%
  data.frame() %>%
  arrange(desc(Freq)) %>%
  rename(gene_symbol = ".")
deg_filtered <- deg_count_df %>%
  filter(Freq == 3)

# filter by bulk DEGs -----------------------------------------------------
## merge bulk rna and protein
deg_bulk_df <- merge(deg_bulk_rna_df, deg_bulk_pro_df, by = c("gene_symbol"), suffixes = c("_rna", "_pro"))
## only keep DEGs with the same direction
deg_bulk_filtered_df <- deg_bulk_df %>%
  filter(!is.na(fdr_rna)) %>%
  filter(!is.na(fdr_pro)) %>%
  filter(fdr_rna < 0.05 & fdr_pro < 0.05) %>%
  filter(number_values >= 20) %>%
  filter((meddiff_exp_rna > 0) == (meddiff_exp_pro > 0)) %>%
  mutate(direction_tumor_vs_normal = ifelse(meddiff_exp_rna > 0, "tumor+", "tumor-"))
deg_filtered <- deg_filtered %>%
  filter(gene_symbol %in% deg_bulk_filtered_df$gene_symbol)
deg_filtered$direction_tumor_vs_normal_bulk <- mapvalues(x = deg_filtered$gene_symbol, from = deg_bulk_filtered_df$gene_symbol, to = as.vector(deg_bulk_filtered_df$direction_tumor_vs_normal))
deg_filtered$direction_tumor_vs_normal_bulk <- as.vector(deg_filtered$direction_tumor_vs_normal_bulk)

# filter pseudotime degs by overlapping with cell type DEGs -----------------------------------
## map the case ids
deg_by_celltype_df$Id_Case <- mapvalues(x = deg_by_celltype_df$aliquot_tumor, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
## filter by case
deg_by_celltype_filtered_df <- deg_by_celltype_df %>%
  filter(Id_Case %in% deg_by_pseudime_df$Id_Run) %>%
  mutate(direction_tumor_vs_normal = ifelse(avg_diff > 0, "tumor+", "tumor-"))
## identify genes always with the same direction
deg_by_celltype_wide_df <- dcast(data = deg_by_celltype_filtered_df, formula = gene_symbol ~ aliquot_tumor, value.var = c("direction_tumor_vs_normal"))
deg_by_celltype_wide_df$direction_summary <- sapply(1:nrow(deg_by_celltype_wide_df), function(i, deg_df) {
  if (any(is.na(deg_by_celltype_wide_df[i,2:5]))) {
    return(NA)
  } else {
    direction_tmp <- unique(unlist(deg_by_celltype_wide_df[i,2:5]))
    if (length(direction_tmp) == 1) {
      return(direction_tmp)
    } else {
      return(NA)
    }
  }
}, deg_df = deg_by_celltype_wide_df)
degs_celltype_same_direction <- deg_by_celltype_wide_df$gene_symbol[!is.na(deg_by_celltype_wide_df$direction_summary)]
deg_filtered <- deg_filtered %>%
  filter(gene_symbol %in% degs_celltype_same_direction)
deg_filtered$direction_tumor_vs_normal_sn <- mapvalues(x = deg_filtered$gene_symbol, from = deg_by_celltype_wide_df$gene_symbol, to = as.vector(deg_by_celltype_wide_df$direction_summary))
deg_filtered$direction_tumor_vs_normal_sn <- as.vector(deg_filtered$direction_tumor_vs_normal_sn)

# filter the direciton of bulk and sn together ----------------------------
deg_filtered <- deg_filtered %>%
  filter(direction_tumor_vs_normal_bulk == direction_tumor_vs_normal_sn)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PT_and_TumorCells_ByCellType.Pseudotim_DE_genes.Combine3Samples.", "Filtered.", run_id, ".tsv")
write.table(x = deg_filtered, file = file2write, quote = F, sep = "\t", row.names = F)

