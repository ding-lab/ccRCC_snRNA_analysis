# Yige Wu @WashU Sep 2019
## merge the human protein atlas data with manually curated marker genes

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# set parameters ----------------------------------------------------------
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input human protein atlas data ------------------------------------------
proteinatlas_tab <- fread(input = "./Ding_Lab/Databases/Human_Protein_Atlas/proteinatlas.tsv", data.table = F)

# input manual curated markers --------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - RCC_Marker_Tab_w.HumanProteinAtlas.20191002.v2.tsv")

# Maybe: to get the kidney-specific genes from human protein atlas ---------------
## at RNA or protein level
## tissue enhanced

# merge the kidney-specific genes -----------------------------------------
genes2print <- unique(gene2cellType_tab$Gene)

# merge human protein atlas with manual curation --------------------------
marker_tab <- data.frame(Gene = genes2print)
marker_tab <- merge(marker_tab, gene2cellType_tab %>%
                      select(Gene, Cell_Type_Abbr, `Structure/Cell Type`, `Sub-Structure/Sub-Type`, Biomarker_Type, Comments, PMIDs, Experimental_Method), by = c("Gene"), all.x = T)
marker_tab <- merge(marker_tab, proteinatlas_tab, by = c("Gene"), all.x = T)
marker_tab <- marker_tab %>%
  arrange(Cell_Type_Abbr)
write.table(x = marker_tab, file = paste0(dir_out, "Kidney_Gene2CellType_Tab_w.HumanProteinAtlas", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp ,".tsv"), quote = F, sep = "\t", row.names = F)
