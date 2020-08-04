# Yige Wu @WashU Aug 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(monocle)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# make colors for tumor segments ------------------------------------------
colors_tumor_segments <- c("T1" = "#FF7F00", "T2" = "#C2A5CF", "T3"  = "#7B3294", "N" = "#008837")

# make color palette for major cell groups -----------------------------------
# cellgroup_colors <- Polychrome::palette36.colors(n = 36)[3:6]
cellgroup_colors <- RColorBrewer::brewer.pal(n = 4, name = "Dark2")[c(4,3,2,1)]
cellgroup_colors <- c(cellgroup_colors, "grey50")
names(cellgroup_colors) <- c("Tumor cells", "Immune",  "Stroma", "Normal epithelial cells", "Unknown")

# make color palette for immune cell groups -----------------------------------
immunecelltype1_colors <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")[c(5,7,8)]
names(immunecelltype1_colors) <- c("Myleoid lineage immune cells", "Lymphoid lineage immune cells",  "Mixed myeloid/lymphoid lineage immune cells")

# make color palette for detailed cell types -----------------------------------
normal_epithelial_colors <- Polychrome::palette36.colors(n = 36)[7:12]
names(normal_epithelial_colors) <- c("Distal convoluted tubule",
                                     "Intercalated cells",
                                     "Principle cells",
                                     "Loop of Henle", 
                                     "Podocytes", 
                                     "Proximal tubule")
# swatch(normal_epithelial_colors)
stroma_colors <- Polychrome::palette36.colors(n = 36)[c(23, 24, 27)]
names(stroma_colors) <- c("Endothelial cells",
                          "Fibroblasts",
                          "Myofibroblasts")
# swatch(stroma_colors)
immune_lymphoid_colors <- colorblind_pal()(8)
# immune_lymphoid_colors <- Polychrome::palette36.colors(n = 36)[c(28, 14:20)]
names(immune_lymphoid_colors) <- c("B-cells", 
                                   "Plasma cells", 
                                   "CD4/CD8 proliferating", 
                                   "CD4+ T-cells",
                                   "Tregs", 
                                   "CD8+ T-cells activated", 
                                   "CD8+ T-cells exhausted",
                                   "NK cells")
# swatch(immune_lymphoid_colors)
immune_myeloid_colors <- Polychrome::palette36.colors(n = 36)[c(13:18)]
names(immune_myeloid_colors) <- c("Basophils", 
                                  "cDC", 
                                  "Macrophages", 
                                  "Macrophages M2b", 
                                  "pDC",
                                  "TRM")
# swatch(immune_myeloid_colors)
immune_mixed_color <- immunecelltype1_colors["Myleoid lineage immune cells"]; names(immune_mixed_color) <- "Mixed myeloid/lymphoid"
immune_colors <- c(immune_lymphoid_colors, immune_myeloid_colors, immune_mixed_color)
tumor_unknown_colors <- cellgroup_colors[c("Tumor cells", "Unknown", "Normal epithelial cells")]
names(tumor_unknown_colors) <- c("Tumor cells", "Unknown", "Normal epithelial cells")
celltype_shorter_colors <- c(tumor_unknown_colors, stroma_colors, normal_epithelial_colors, immune_colors)

# write ouput -------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_color_palettes.RData")
save(celltype_shorter_colors, cellgroup_colors, colors_tumor_segments, file = file2write)

