# Yige Wu @WashU March 2020
## make barcode to cell type mapping table for the integrated dataset

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode-cell-type annotation
barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/other/map_celltype_to_barcode/20200320.v1/30_aliquot_integration.barcode2celltype.20200320.v1.tsv", data.table = F)
## set group "patient" cells to extract
snRNA_cellgroup2process <- c("Endothelial cells", "B-cells", "CD4\\+ T-cells", "CD8\\+ T-cells", "NK cells", "Fibroblasts", "Macrophages", "DC", "Nephron_Epithelium")
xCell_parent_cells2process <- c("Endothelial cells", "B-cells", "CD4+ T-cells", "CD8+ T-cells", "NK cells", "Fibroblasts", "Macrophages", "DC", "Epithelial cells")

# count by cell type group ------------------------------------------------
cellgroup_count_df <- NULL
# cellgroup_tmp <- "NK cells"
for (cellgroup_tmp in snRNA_cellgroup2process) {
  barcode2celltype_tmp_df <- barcode2celltype_df %>%
    filter(grepl(pattern = cellgroup_tmp, x = Most_Enriched_Cell_Group, ignore.case = F) |grepl(pattern = cellgroup_tmp, x = Most_Enriched_Cell_Type1, ignore.case = F) | grepl(pattern = cellgroup_tmp, x = Most_Enriched_Cell_Type2, ignore.case = F)| grepl(pattern = cellgroup_tmp, x = Most_Enriched_Cell_Type3, ignore.case = F)| grepl(pattern = cellgroup_tmp, x = Most_Enriched_Cell_Type4, ignore.case = F))
  cellgroup_count_tmp_df <- barcode2celltype_tmp_df %>%
    select(orig.ident) %>%
    table() %>%
    as.data.frame() %>%
    mutate(cellgroup = cellgroup_tmp)
  cellgroup_count_df <- rbind(cellgroup_count_tmp_df, cellgroup_count_df)
}

# add additional xCell for merging in the future --------------------------
## rename columns
cellgroup_stat_df <- cellgroup_count_df %>%
  rename(Aliquot.snRNA = '.') %>%
  rename(CellgroupFreq.snRNA = Freq) %>%
  rename(Cellgroup.snRNA = cellgroup)
cellgroup_stat_df$CellgroupFreq.snRNA <- as.vector(cellgroup_stat_df$CellgroupFreq.snRNA)
## calculate total cell count per sample
total_count_df <- barcode2celltype_df %>%
  select(orig.ident) %>%
  table() %>%
  as.data.frame() %>%
  rename(Aliquot.snRNA = '.') %>%
  rename(AliquotTotalFreq.snRNA = Freq)
total_count_df$AliquotTotalFreq.snRNA <- as.vector(total_count_df$AliquotTotalFreq.snRNA)
## add total cell count per sample
cellgroup_stat_df <- merge(cellgroup_stat_df, total_count_df, by = c("Aliquot.snRNA"), all.x = T)
cellgroup_stat_df$CellgroupFrac.snRNA <- cellgroup_stat_df$CellgroupFreq.snRNA/cellgroup_stat_df$AliquotTotalFreq.snRNA
## add xcell info
cellgroup_stat_df$Cellgroup.xCell <- mapvalues(x = cellgroup_stat_df$Cellgroup.snRNA, from = snRNA_cellgroup2process, to = xCell_parent_cells2process)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "snRNA_cellgroup_stats_for_xCell_comparison.", run_id, ".tsv")
write.table(x = cellgroup_stat_df, file = file2write, quote = F, sep = "\t", row.names = F)
