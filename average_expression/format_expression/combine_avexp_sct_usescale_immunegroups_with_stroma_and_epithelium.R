# Yige Wu @WashU Aug 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
avgexp_immune_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/averageexpression_immune_sct_usescale_byimmunegroup_on_katmai/20200908.v1/averageexpression.SCT.byCell_group.immune.31_aliquot_integration.20200908.v1.tsv")
avgexp_bycelltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_bycelltypeshorter_on_katmai/20200904.v1/averageexpression_SCT_bycelltype.shorter.31_aliquot_integration.20200904.v1.tsv")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_immune_cellgroups/20200908.v1/31Aliquot.Barcode2CellType.Immune_Grouped.20200908.v1.tsv", data.table = F)

# get the non-immune cell types to keep ------------
## count cells
celltype_keep_df <- barcode2celltype_df %>%
  select(Cell_type.shorter) %>%
  table() %>%
  as.data.frame() %>%
  rename(Cell_type.shorter = '.')
## get celltype-to-cellgroup map
celltype2cellgroup_df <- barcode2celltype_df %>%
  select(Cell_type.shorter, Cell_group.shorter, Cell_group.shorter) %>%
  unique() %>%
  mutate(colname_celltype = gsub(x = Cell_type.shorter, pattern = " |\\+|\\/|\\-", replacement = "."))
## map cell group
celltype_keep_df$Cell_group.shorter <- mapvalues(x = celltype_keep_df$Cell_type.shorter, from = celltype2cellgroup_df$Cell_type.shorter, to = as.vector(celltype2cellgroup_df$Cell_group.shorter))
## get cell type names to keep
celltype_keep_df <- celltype_keep_df %>%
  mutate(Keep = (Freq >= 50 & Cell_group.shorter != "Immune" & !(Cell_type.shorter %in% c("Unknown", "Tumor-like cells")))) %>%
  mutate(colname_celltype = gsub(x = Cell_type.shorter, pattern = " |\\+|\\/|\\-", replacement = "."))
## subset average expression
### filtr the rows
avgexp_nonimmune_df <- avgexp_bycelltype_df %>%
  rename(gene = V1)
### remove teh prefix from the column names
data_col_names <- colnames(avgexp_nonimmune_df)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
data_col_names.changed
### rename the data frame
colnames(avgexp_nonimmune_df) <- c("gene", data_col_names.changed)
### subset
avgexp_nonimmune_df <- avgexp_nonimmune_df[, c("gene", celltype_keep_df$colname_celltype[celltype_keep_df$Keep])]

# get the immune cell types to keep ------------
## get immune cell types to keep
immunecelltype_keep_df <- barcode2celltype_df %>%
  select(Cell_group.immune) %>%
  unique() %>%
  mutate(Keep = ifelse(Cell_group.immune %in% c("", "others"), F, T)) %>%
  mutate(colname_celltype = gsub(x = Cell_group.immune, pattern = " |\\+|\\/|\\-", replacement = "."))
## subset average expression
### filtr the rows
avgexp_immune_df <- avgexp_immune_df %>%
  rename(gene = V1)
### remove teh prefix from the column names
data_col_names <- colnames(avgexp_immune_df)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
data_col_names.changed
### rename the data frame
colnames(avgexp_immune_df) <- c("gene", data_col_names.changed)
### subset
avgexp_immune_df <- avgexp_immune_df[, c("gene", immunecelltype_keep_df$colname_celltype[immunecelltype_keep_df$Keep])]

# combine data frames -----------------------------------------------------
head(avgexp_immune_df$gene)
head(avgexp_nonimmune_df$gene)
avgexp_df <- cbind(avgexp_immune_df, avgexp_nonimmune_df[,-1])
## combine the cell type to cell group table
celltype2cellgroup_keep_df <- rbind(celltype_keep_df %>%
                                      filter(Keep) %>%
                                      rename(Cell_type = Cell_type.shorter) %>%
                                      select(Cell_type, colname_celltype, Cell_group.shorter),
                                    immunecelltype_keep_df %>%
                                      filter(Keep) %>%
                                      rename(Cell_type = Cell_group.immune) %>%
                                      mutate(Cell_group.shorter = "Immune") %>%
                                      select(Cell_type, colname_celltype, Cell_group.shorter))
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Combined.", "Average_Expression.", "SCT.", "UseScale.", run_id, ".tsv")
write.table(x = avgexp_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Combined.", "CellTypes2CellGroup.", run_id, ".tsv")
write.table(x = celltype2cellgroup_keep_df, file = file2write, quote = F, sep = "\t", row.names = F)


