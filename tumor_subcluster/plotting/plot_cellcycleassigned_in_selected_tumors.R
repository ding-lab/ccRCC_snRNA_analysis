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
packages = c(
  "RCurl"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input seurat paths
paths_srat_df <- fread(data.table = F, input = "./Data_Freezes/V1/snRNA/Tumor_Cell_Reclustered/Paths_TumorCellOnlyReclustered_SeuratObject.20201127.v1.tsv")
## input the barcode-manualsubcluster info
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20200902.v1/scrublet.run20200902_adj_cutoff.united_outputs.tsv", data.table = F)
## retreive human cell cycle genes
cell_cycle_genes <- fread(input = "../ccRCC_Drug/Resources/Analysis_Results/dependencies/make_cellcycle_human_mouse_genes/20200408.v1/cell_cycle_human_mouse_genes.20200408.v1.tsv", data.table = F)

# preprocess --------------------------------------------------------------
cell_cycle_genes <- cell_cycle_genes %>%
  filter(species == "human")
g2m_feature_names <- cell_cycle_genes$gene_name[cell_cycle_genes$phase == "G2/M"]
s_feature_names <- cell_cycle_genes$gene_name[cell_cycle_genes$phase == "S"]

# process each sample -----------------------------------------------------
easyid_tmp <- "C3L-00583-T1"
easyid_tmp <- "C3L-01313-T1"
phase_df <- NULL
for (easyid_tmp in unique(paths_srat_df$Aliquot.snRNA.WU)) {
  ## input seurat object,
  path_srat_katmai <- paths_srat_df$Path_katmai[paths_srat_df$Aliquot.snRNA.WU == easyid_tmp]
  path_srat_relative <- gsub(x = path_srat_katmai, pattern = "\\/diskmnt\\/Projects\\/ccRCC_scratch\\/ccRCC_snRNA\\/", replacement = "./")
  srat <- readRDS(file = path_srat_relative)
  DefaultAssay(srat) <- "RNA"
  ## subset
  scrublets_df <- barcode2scrublet_df %>%
    filter(Aliquot_WU == easyid_tmp) %>%
    filter(predicted_doublet)
  barcodes_keep <- rownames(srat@meta.data)
  barcodes_keep <- barcodes_keep[!(barcodes_keep %in% scrublets_df$Barcodes)]
  srat <- subset(x = srat, cells = barcodes_keep)
  ## assign cell cycle scores
  srat <- CellCycleScoring(srat, s.features = s_feature_names, g2m.features = g2m_feature_names, set.ident = TRUE)
  ## score phase info
  phase_tmp_df <- srat@meta.data
  phase_tmp_df$barcode <- rownames(phase_tmp_df)
  phase_df <- rbind(phase_df, phase_tmp_df)
  ## plot
  p <- DimPlot(srat, group.by= "Phase")
  ## write output
  file2write <- paste0(dir_out, easyid_tmp, ".png")
  png(file2write, width = 1000, height = 1000, res = 150)
  print(p)
  dev.off()
}

# write phase info --------------------------------------------------------
file2write <- paste0(dir_out, "TumorCellPhase.", "DoubletRemoved", ".tsv")
write.table(x = phase_df, file = file2write, sep = "\t", row.names = F, quote = F)


