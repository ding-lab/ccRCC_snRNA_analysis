# Yige Wu @WashU Sep 2020

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
## input the integrated data
path_rds <- "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/subset_C3N-01200_tumorlikecells_and_recluster/20200910.v1/TumorLikeCells.Reclustered.20200910.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-manualsubcluster info
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200904.v1/31Aliquot.Barcode2CellType.20200904.v1.tsv", data.table = F)

# add cell type to the Seurat meta data---------------------------------------------
BC <- srat@meta.data %>% rownames
## get original barcode
srat@meta.data$original_barcode <- BC %>% strsplit("_") %>% lapply("[[",1) %>% unlist
head(srat@meta.data$original_barcode)
## make combined id for the seurat meta data
srat@meta.data$id_aliquot_barcode <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$original_barcode)
head(srat@meta.data$id_aliquot_barcode)
## make combined id for the barcode2celltype table
barcode2celltype_df <- barcode2celltype_df %>%
  filter(orig.ident %in% unique(srat@meta.data$orig.ident))
barcode2celltype_df$id_aliquot_barcode <- paste0(barcode2celltype_df$orig.ident, "_", barcode2celltype_df$individual_barcode)
## map cell type shorter
srat@meta.data$Cell_type.shorter <- mapvalues(x = srat@meta.data$id_aliquot_barcode, from = barcode2celltype_df$id_aliquot_barcode, to = as.vector(barcode2celltype_df$Cell_type.shorter))

# Dimplot -----------------------------------------------------------------
p <- DimPlot(object = srat, group.by = "Cell_type.shorter")
p
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Celltype.", "C3N-01200", ".png")
png(file2write, width = 900, height = 800, res = 150)
print(p)
dev.off()

