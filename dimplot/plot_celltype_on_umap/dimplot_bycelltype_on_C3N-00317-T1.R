# Yige Wu @WashU Mar 2021

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
srat <- readRDS(file = "./Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0012280004/pf3000_fmin200_fmax10000_cmin3000_cmax10000_mito_max0.1/CPT0012280004_processed.rds")
print("Finish reading RDS file")
## input the barcode-manualsubcluster info
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/individual_sample/map_celltype_for_C3N-00317-T1/20210305.v1/C3N-00317-T1.Barcode2CellType.20210305.v1.tsv", data.table = F)

# add cell type to the Seurat meta data---------------------------------------------
BC <- srat@meta.data %>% rownames
## map cell type shorter
srat@meta.data$Cell_group13 <- mapvalues(x = BC, from = barcode2celltype_df$individual_barcode, to = as.vector(barcode2celltype_df$Cell_group13))

# Dimplot -----------------------------------------------------------------
p <- DimPlot(object = srat, group.by = "Cell_group13")
p <- p + scale_color_manual(values = colors_cellgroup13)
p
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cell_group13",".png")
png(file2write, width = 900, height = 800, res = 150)
print(p)
dev.off()

