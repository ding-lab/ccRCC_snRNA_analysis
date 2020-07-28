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
library(ggplot2)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Resources/Analysis_Results/integration/31_aliquot_integration/31_aliquot_integration_without_anchoring/20200727.v1/31_aliquot_integration_without_anchoring.20200727.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
path_barcode2celltype <- "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells_with_patch/20200720.v1/31AliquotIntegration.Barcode2CellType.TumorManualCluster.20200720.v1.tsv"
barcode2celltype_df <- fread(input = path_barcode2celltype, data.table = F)
## input DEG for each cell type
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200406.v1.tsv")
## specify minimal percentage of cell in any clustering expressing the genes to show
min.exp.pct <- 20
min.exp.pct.malignantmarkers <- 10

# modify srat object ------------------------------------------------------
## remove unknown cell group 
metadata_df <- srat@meta.data
metadata_df$integrated31_barcode <- rownames(metadata_df)
metadata_df <- metadata_df %>%
  mutate(individual_barcode = str_split_fixed(string = integrated31_barcode, pattern = "_", n = 3)[,1])
barcode2celltype_df <- merge(barcode2celltype_df, metadata_df, by = c("orig.ident", "individual_barcode"), all.x = T)
barcode2celltype_filtered_df <- barcode2celltype_df %>%
  filter(Most_Enriched_Cell_Group == "Nephron_Epithelium")
barcode2celltype_filtered_df$Cell_type.detailed %>% table()
barcodes2process <- barcode2celltype_filtered_df$integrated31_barcode
srat <- subset(srat, cells = barcodes2process)
## change meta data
srat@meta.data$Cell_type.detailed <- mapvalues(x = rownames(srat@meta.data), from = barcode2celltype_df$integrated31_barcode, to = as.vector(barcode2celltype_df$Cell_type.detailed))
srat@meta.data$Cell_type.detailed %>% table()
## change ident
Idents(srat) <- "Cell_type.detailed"

# get genes to plot -------------------------------------------------------
## make the data frame for the dotplot
## get the genes within the cell type marker table
genes2plot <-  intersect(gene2celltype_df$Gene[gene2celltype_df$Cell_Type_Group %in% c("Nephron_Epithelium", "Malignant_Nephron_Epithelium")], srat@assays$RNA@counts@Dimnames[[1]])
genes2plot <- unique(genes2plot)
## get the pct expressed for each gene in each cluster
p <- DotPlot(object = srat, features = genes2plot, col.min = 0)
plot_data <- p$data
## transform the dataframe to matrix to better filter out genes with too low expressin
plot_matrix <- dcast(data = plot_data, formula = features.plot ~ id, value.var = "pct.exp")
plot_matrix %>% head()
## replot with the filtered genes plus malignant cell marker genes
### filter for genes that are expressed in > min.exp.pct % of one cluster at least
genes2plot_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > min.exp.pct) >= 1, "features.plot"])
malignant_markers <- as.vector(gene2celltype_df$Gene[gene2celltype_df$Cell_Type_Group == "Malignant_Nephron_Epithelium"])
malignant_markers_filtered <- as.vector(plot_matrix[(plot_matrix$features.plot %in% malignant_markers) & (rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > min.exp.pct.malignantmarkers) >= 1), "features.plot"])
genes2plot_filtered <- unique(c(genes2plot_filtered, malignant_markers_filtered))

# plot --------------------------------------------------------------------
## make the dotplot
cat("###########################################\n")
cat("Dotplot now\n")
p <- DotPlot(object = srat, features = genes2plot_filtered, col.min = 0)
p$data$id <- factor(x = p$data$id, levels = rev(c("Tumor cells", "Proximal tubule", "Loop of Henle", "Distal convoluted tubule", "Principle cells", "Intercalated cells", "Podocytes")))
# p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
# p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
# p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
# p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
# p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
# p <- p  + RotatedAxis()
# p <- p + facet_grid(.~gene_cell_type1 + gene_cell_type2 + gene_cell_type3, scales = "free", space = "free", drop = T)
# p <- p + theme(panel.spacing = unit(0, "lines"),
#                strip.background = element_blank(),
#                panel.border = element_rect(colour = "black"),
#                strip.text.x = element_text(angle = 90, vjust = 0.5),
#                strip.placement = "outside")
# p <- p + theme(panel.grid.major = element_line(colour = "grey50"))
p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold"),
               axis.text.y = element_text(face = "bold"))
p <- p + theme(axis.title.x = element_blank(),
               axis.title.y = element_blank())
cat("Finished Dotplot\n")
cat("###########################################\n")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Dotplot_CellTypeMarkers_Exp", ".pdf")
pdf(file = file2write, width = 20, height = 5, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, "Dotplot_CellTypeMarkers_Exp", ".RDS")
saveRDS(object = p, file = file2write, compress = T)
