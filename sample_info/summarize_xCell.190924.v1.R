# Yige Wu @WashU Sep 2019
## for summarizing xCell estimate of immune cell + stroma cell fractions from CPTAC ccRCC manuscript

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# Input xCell stats -------------------------------------------------------
xcell_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S7.xlsx", sheet = "xCell Signatures", skip = 2)
## because there is a row in character format so the rest of the data frame is in character format, need to change it double format
xcell_value_tab <- xcell_tab %>%
  filter(Samples != "Immune Group")

xcell_value_melt_tab <- melt(xcell_value_tab, id.vars = c("Samples"))
colnames(xcell_value_melt_tab) <- c("Cells", "SampID", "Value")
xcell_value_melt_tab$Value <- as.numeric(as.vector(xcell_value_melt_tab$Value))
xcell_value_melt_tab %>% head()

# input table of cell family from xCell -----------------------------------
xcell_cell_family_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/xCell/Dev_scripts/cells_families.txt", data.table = F)
xcell_cell_family_tab %>% head()
xcell_cell_family_tab %>%
  select(Group) %>%
  table()
xcell_cell_family_tab %>%
  select(Family) %>%
  table()
xcell_cell_family_tab %>%
  select(Type) %>%
  table()
xcell_cell_family_tab %>%
  select(Family2) %>%
  table()


# check which xCell groups to summarize -----------------------------------
xcell_cell_family_tab %>%
  select(Cells) %>%
  table()
xcell_cell_family_tab %>%
  select(Group) %>%
  table()

xcell_cells2sum <- list()
xcell_cells2sum[["MacroMono"]] <- xcell_cell_family_tab %>%
  filter(Group %in% c("Macrophages", "Monocytes", "Macrophages M1", "Macrophages M2") | Cells %in% c("Macrophages", "Monocytes", "Macrophages M1", "Macrophages M2")) %>%
  select(Cells) %>%
  unlist(use.names = F)

xcell_cells2sum[["Tcell"]] <- xcell_cell_family_tab %>%
  filter(Group %in% c("CD4+ memory T-cells", "CD4+ T-cells", "CD8+ T-cells") | Cells %in% c("CD4+ memory T-cells", "CD4+ T-cells", "CD8+ T-cells")) %>%
  select(Cells) %>%
  unlist(use.names = F)

xcell_cells2sum[["pDC"]] <- xcell_cell_family_tab %>%
  filter(Group %in% c("pDC") | Cells %in% c("pDC")) %>%
  select(Cells) %>%
  unlist(use.names = F)


xcell_cells2sum[["Endothelial"]] <- xcell_cell_family_tab %>%
  filter(Group %in% c("Endothelial cells") | Cells %in% c("Endothelial cells")) %>%
  select(Cells) %>%
  unlist(use.names = F)

xcell_cells2sum[["Epithelial"]] <- xcell_cell_family_tab %>%
  filter(Group %in% c("Epithelial cells") | Cells %in% c("Epithelial cells")) %>%
  select(Cells) %>%
  unlist(use.names = F)


xcell_cells2sum

# merge value with cell family info ---------------------------------------
xcell_value_melt_tab <- merge(xcell_value_melt_tab, xcell_cell_family_tab, by = c("Cells"), all.x = T)
xcell_value_melt_tab %>% head()

# summarize and write table -----------------------------------------------
class(xcell_value_melt_tab$Value)

xcell_value_sum_tab <- data.frame(SampID = unique(xcell_value_melt_tab$SampID))
for (cell_cluster_tmp in names(xcell_cells2sum)) {
  xcell_value_tmp_df <- xcell_value_melt_tab %>%
    filter(Cells %in% xcell_cells2sum[[cell_cluster_tmp]]) %>%
    group_by(SampID) %>%
    summarize(xcell_value_tmp = sum(Value))
  colnames(xcell_value_tmp_df)[2] <- cell_cluster_tmp
  xcell_value_sum_tab <- merge(xcell_value_sum_tab, xcell_value_tmp_df, by = c("SampID"), all.x = T)
}


# write -------------------------------------------------------------------
version_tmp <- 1
write.table(x = xcell_value_sum_tab, file = paste0(makeOutDir(), "xCell_value_sum", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp , ".tsv"), quote = F, sep = "\t", row.names = F)
