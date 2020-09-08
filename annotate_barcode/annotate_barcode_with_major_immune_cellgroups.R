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
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200904.v1/31Aliquot.Barcode2CellType.20200904.v1.tsv", data.table = F)

# group detailed immune cell types into major immune cell groups ----------
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(Cell_group.immune = ifelse(Cell_group.shorter == "Immune",
                                    ifelse(Cell_type1 == "Myleoid lineage immune cells",
                                           ifelse(Cell_type3 == "Macrophages", 
                                                  "Macrophages",
                                                  ifelse(Cell_type3 == "DC", 
                                                         "DC",
                                                         "others")),
                                           ifelse(Cell_type2 == "NK cells",
                                                  "NK cells",
                                                  ifelse(Cell_type2 == "T-cells",
                                                         ifelse(Cell_type3 == "CD4+ T-cells",
                                                               "CD4+ T-cells",
                                                               ifelse(Cell_type3 == "CD8+ T-cells", 
                                                                      "CD8+ T-cells", 
                                                                      "others")),
                                                         ifelse(Cell_type2 == "B-cells", 
                                                                "B-cells", 
                                                                "others")))),
                                    ""))

barcode2celltype_df %>%
  filter(Cell_group.shorter == "Immune") %>%
  select(Cell_type.detailed, Cell_group.immune) %>%
  unique() %>%
  arrange(Cell_group.immune)

# write output ------------------------------------------------------------
## final check up
nrow(barcode2celltype_df) # 138547
## 94846 tumor cells, 989 unknown cells, 213 tumor-like epithelial cells
file2write <- paste0(dir_out, "31Aliquot.Barcode2CellType.Immune_Grouped.", run_id, ".tsv")
write.table(x = barcode2celltype_df, file = file2write, quote = F, sep = "\t", row.names = F)





