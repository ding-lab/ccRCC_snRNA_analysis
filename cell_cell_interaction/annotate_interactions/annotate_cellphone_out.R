# Yige Wu @WashU Sep 2020
## reference of the cellphonedb output: https://www.cellphonedb.org/documentation

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
## input cellphonedb output
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/filter_interactions/filter_cellphonedb_out/20200925.v1/cell.phone.res.total.run20200818.filtered.txt")
## input cell type annotation
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20200917.v2/31Aliquot.Barcode2CellType.20200917.v2.tsv")

# add sample type and rank of the value within the same case --------------
cellphone_df <- cellphone_df %>%
  mutate(Sample_type = ifelse(grepl(x = Easy_id, pattern = "-N"), "Normal", "Tumor")) %>%
  group_by(Easy_id) %>%
  mutate(rank_sig_mean_same_case = order(order(value, decreasing = T)))

# edit gene_a and gene_b --------------------------------------------------
cellphone_df$gene_a <- sapply(X = cellphone_df$interacting_pair, FUN = function(x) {
  text_vec <- str_split(string = x, pattern = "_")[[1]]
  if (length(text_vec) == 2) {
    gene_a <- text_vec[1]
    return(gene_a)
  } else {
    gene_a <- paste0(text_vec[-length(text_vec)], collapse = "_")
    return(gene_a)
  }
})
cellphone_df$gene_b <- sapply(X = cellphone_df$interacting_pair, FUN = function(x) {
  text_vec <- str_split(string = x, pattern = "_")[[1]]
  if (length(text_vec) == 2) {
    gene_b <- text_vec[2]
    return(gene_b)
  } else {
    gene_b <- text_vec[length(text_vec)]
    return(gene_b)
  }
})

# annotate receptor gene and ligand gene ----------------------------------
cellphone_df$receptor_b[cellphone_df$is_integrin] <- T
cellphone_new_df <- cellphone_df %>%
  mutate(is_ligand_receptor = ((receptor_a & !receptor_b) | (receptor_b & !receptor_a))) %>%
  mutate(gene.source = ifelse(receptor_a, gene_b, gene_a)) %>%
  mutate(gene.target = ifelse(receptor_a, gene_a, gene_b)) %>%
  mutate(receptor.source = ifelse(gene.source == gene_a, receptor_a, receptor_b)) %>%
  mutate(receptor.target = ifelse(gene.target == gene_a, receptor_a, receptor_b))

# correct source gene and target gene -------------------------------------
# ERBB3_NRG1
## ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4712818/
## ref: Here, we demonstrate neuregulin 1 (NRG1) expression and secretion by human decidual stromal cells. Stimulation of human primary trophoblasts with exogenous NRG1 induced phosphorylation of ErbB2, ErbB3 and related downstream effectors. Co-immunoprecipitation experiments confirmed the formation of ErbB2–ErbB3 dimers upon ligand engagement. Along this line, receptor knockdown and ErbB3 neutralization strongly diminished NRG1-dependent activation of the signalling complex.
cellphone_new_df$gene.source[cellphone_new_df$interacting_pair == "ERBB3_NRG1"] <- "NRG1"
cellphone_new_df$gene.target[cellphone_new_df$interacting_pair == "ERBB3_NRG1"] <- "ERBB3"
cellphone_new_df$is_ligand_receptor[cellphone_new_df$interacting_pair == "ERBB3_NRG1"] <- T
cellphone_new_df$gene.source[cellphone_new_df$interacting_pair == "NRG1_LGR4"] <- "NRG1"
cellphone_new_df$gene.target[cellphone_new_df$interacting_pair == "NRG1_LGR4"] <- "LGR4"
cellphone_new_df$is_ligand_receptor[cellphone_new_df$interacting_pair == "NRG1_LGR4"] <- T
# EGFR_HBEGF
## ref: https://www.genecards.org/cgi-bin/carddisp.pl?gene=HBEGF
cellphone_new_df$gene.source[cellphone_new_df$interacting_pair == "EGFR_HBEGF"] <- "HBEGF"
cellphone_new_df$gene.target[cellphone_new_df$interacting_pair == "EGFR_HBEGF"] <- "EGFR"
cellphone_new_df$receptor.source[cellphone_new_df$interacting_pair == "EGFR_HBEGF"] <- F
cellphone_new_df$is_ligand_receptor[cellphone_new_df$interacting_pair == "EGFR_HBEGF"] <- T
cellphone_new_df <- cellphone_new_df %>%
  mutate(interaction.source2target = ifelse(is_ligand_receptor, 
                                            paste0(gene.source, "->", gene.target),
                                            paste0(gene.source, "&", gene.target)))


# annotate cell type source and target ------------------------------------
cellphone_new_df <- cellphone_new_df %>%
  mutate(Cell_type.source = ifelse(gene.source == gene_a, Cell_type1, Cell_type2)) %>%
  mutate(Cell_type.target = ifelse(gene.source == gene_a, Cell_type2, Cell_type1)) %>%
  mutate(celltypes.source2target = ifelse(is_ligand_receptor, 
                                            paste0(Cell_type.source, "->", Cell_type.target),
                                            paste0(Cell_type.source, "&", Cell_type.target)))

# annotate cell type to cell groups ---------------------------------------
celltype2cellgroup_df <- barcode2celltype_df %>%
  select(Cell_type.shorter, Cell_group7) %>%
  unique()
celltype2cellgroup_df <- rbind(celltype2cellgroup_df,
                               data.frame(Cell_type.shorter = c("CD4 CTL", "CD4 T-cells", "CD4 T-cells activated", "CD4 T-cells naive", "CD8 CTL", 
                                                                "CD8 CTL exhausted", "CD8 T-cells preexhausted", "Macrophages", "Macrophages proliferating", 
                                                                "NK cells strong", "NK cells weak", "TRM", "Plasma"),
                                          Cell_group7 = "Immune"))
cellphone_new_df$Cell_group.source <- mapvalues(x = cellphone_new_df$Cell_type.source, from = celltype2cellgroup_df$Cell_type.shorter, to = as.vector(celltype2cellgroup_df$Cell_group7))
# table(cellphone_df$Cell_group1)
cellphone_new_df$Cell_group.target <- mapvalues(x = cellphone_new_df$Cell_type.target, from = celltype2cellgroup_df$Cell_type.shorter, to = as.vector(celltype2cellgroup_df$Cell_group7))

# group paired cell types -------------------------------------------------
table(cellphone_new_df$Cell_group.source)
cellphone_new_df <- cellphone_new_df %>%
  mutate(paired_cellgroups.detailed = ifelse(is_ligand_receptor, 
                                             paste0(ifelse(Cell_group.source == "Tumor cells", "Tumor", Cell_group.source),"->", ifelse(Cell_group.target == "Tumor cells", "Tumor", Cell_group.target)),
                                             paste0(ifelse(Cell_group.source == "Tumor cells", "Tumor", Cell_group.source),"&", ifelse(Cell_group.target == "Tumor cells", "Tumor", Cell_group.target))))
table(cellphone_new_df$paired_cellgroups.detailed)
paired_cellgroups.general_vec <- vector(mode = "character", length = nrow(cellphone_new_df))
paired_cellgroups.general_vec[cellphone_new_df$paired_cellgroups.detailed %in% c("Stroma->Tumor", "Tumor->Stroma", "Stroma&Tumor", "Tumor&Stroma")] <- "Tumor&Stroma"
paired_cellgroups.general_vec[cellphone_new_df$paired_cellgroups.detailed %in% c("Immune->Tumor", "Tumor->Immune", "Immune&Tumor", "Tumor&Immune")] <- "Tumor&Immune"
paired_cellgroups.general_vec[cellphone_new_df$paired_cellgroups.detailed %in% c("Immune->Stroma", "Stroma->Immune", "Immune&Stroma", "Stroma&Immune")] <- "Stroma&Immune"
paired_cellgroups.general_vec[cellphone_new_df$paired_cellgroups.detailed %in% c("Tumor->Tumor", "Tumor&Tumor")] <- "Tumor&Tumor"
paired_cellgroups.general_vec[cellphone_new_df$paired_cellgroups.detailed %in% c("Stroma->Stroma", "Stroma&Stroma")] <- "Stroma&Stroma"
paired_cellgroups.general_vec[cellphone_new_df$paired_cellgroups.detailed %in% c("Immune->Immune", "Immune&Immune")] <- "Immune&Immune"
paired_cellgroups.general_vec[paired_cellgroups.general_vec == ""] <- "Others"
cellphone_new_df$paired_cellgroups.general <- paired_cellgroups.general_vec
table(cellphone_new_df$paired_cellgroups.general)

# add rank by the cellgroup pairs -----------------------------------------
cellphone_new_df <- cellphone_new_df %>%
  group_by(Easy_id, paired_cellgroups.general) %>%
  mutate(rank_sig_mean.paired_cellgroups.general = order(order(value, decreasing = T)))

# remove old columns ------------------------------------------------------
cellphone_new_df <- cellphone_new_df %>%
  select(-receptor_a) %>%
  select(-receptor_b) %>%
  select(-Cell_type1) %>%
  select(-Cell_type2) %>%
  select(-gene_a) %>%
  select(-gene_b)
 

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cell.phone.res.total.run20200818.filtered.formatted.txt")
write.table(x = cellphone_new_df, file = file2write, quote = F, sep = "\t", row.names = F)


