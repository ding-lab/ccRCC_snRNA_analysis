# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Resources/Analysis_Results/integration/30_aliquot_integration/docker_run_integration/20200212.v3/30_aliquot_integration.20200212.v3.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
path_barcode2celltype <- "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells/20200713.v1/30AliquotIntegration.Barcode2CellType.TumorManualCluster.20200713.v1.tsv"
barcode2celltype_df <- fread(input = path_barcode2celltype, data.table = F)
table(barcode2celltype_df$Cell_type.shorter)
## input DEG for each cell type
deg_df <- fread(input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/findallmarker_wilcox_cellgroup_on_katmai/20200714.v1/findallmarkers_roc_bycellgroup.pos.20200714.v1.tsv", data.table = F)
## input hif targets
hif_targets_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_hif_targets/20200428.v1/HIF_Target_Genes.20200428.v1.tsv")

# modify srat object ------------------------------------------------------
## get barcodes to process
barcode2celltype_filtered_df <- barcode2celltype_df %>%
  filter(Cell_group != "Unknown")
## remove unknown cell group 
barcodes2process <- barcode2celltype_filtered_df$integrated_barcode
srat <- subset(srat, cells = barcodes2process)
## change meta data
srat@meta.data$Cell_type.shorter <- mapvalues(x = rownames(srat@meta.data), from = barcode2celltype_df$integrated_barcode, to = as.vector(barcode2celltype_df$Cell_type.shorter))
## change ident
Idents(srat) <- "Cell_type.shorter"

# prepare plotting parameters ---------------------------------------------
## get genes to plot
genes2plot <- unique(hif_targets_df$target_genesymbol)
genes2plot <- intersect(genes2plot, unique(deg_df$gene))
genes2plot
## celltype to cell type category
celltype_cat_df <- barcode2celltype_filtered_df %>%
  select(Cell_type.shorter, Cell_group) %>%
  unique()

# plot --------------------------------------------------------------------
## actual plotting
DefaultAssay(object2plot) <- "RNA"
p <- DotPlot(object = srat, features = genes2plot, col.min = 0)
p$data$Cell_group <- mapvalues(x = p$data$id, from = celltype_cat_df$Cell_type.shorter, to = celltype_cat_df$Cell_group)
p$data$gene_celltype_exp_cat <- mapvalues(x = p$data$features.plot, from = genes2plot, to = gene_celltype_exp_cats)
p$data$gene_celltype_exp_cat <- factor(p$data$gene_celltype_exp_cat, levels = c("Malignant_Expressed", "Stromal_Expressed", "Immune_Expressed"))
p <- p  + RotatedAxis()
p <- p + facet_grid(Cell_group + id~gene_celltype_exp_cat, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background.y = element_rect(colour = "black", fill = "white"),
               strip.background.x = element_rect(colour = "black", fill = "white"),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 12, face = "bold"),
               strip.placement = "outside")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Dotplot_HIF_Downstream_Exp.", ".png")
png(file = file2write, width = 2500, height = 1300, res = 150)
print(p)
dev.off()


