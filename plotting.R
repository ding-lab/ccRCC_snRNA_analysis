# especially for plotting #

# source ------------------------------------------------------------------
packages = c(
  "ggplot2",
  "ggrepel",
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer",
  "ggthemes",
  "rcartocolor",
  "Polychrome"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
cartocolors_df <- rcartocolor::cartocolors


# make colors for tumor segments ------------------------------------------
colors_tumor_segments <- c("T1" = "#FF7F00", "T2" = "#C2A5CF", "T3"  = "#7B3294", "N" = "#008837")

# make color palette for major cell groups -----------------------------------
# cellgroup_colors <- Polychrome::palette36.colors(n = 36)[3:6]
cellgroup_colors <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")[c(4, 2, 6, 3, 1)]
cellgroup_colors <- c(cellgroup_colors, "grey50")
names(cellgroup_colors) <- c("Tumor cells", "Tumor-like epithelial cells", "Normal epithelial cells", "Immune",  "Stroma",  "Unknown")
# swatch(cellgroup_colors)

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
# save(cellgroup_colors, celltype_shorter_colors, file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Figures/r_colorpalette.RData")

# make color palette for variant class ------------------------------------
# rcartocolor::display_carto_all()
cartocolors_temps <- cartocolors_df[cartocolors_df$Name == "Temps", "n7"][[1]]
cartocolors_tropic <- cartocolors_df[cartocolors_df$Name == "Tropic", "n7"][[1]]
variant_class_colors <- c(cartocolors_temps[1:4], 
                          cartocolors_tropic[4], 
                          cartocolors_temps[c(5,7)],
                          "white")
names(variant_class_colors) <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", 'Splice_Site', 
                                 "Silent", 
                                 "Missense_Mutation", "In_Frame_Ins",
                                 "None")


# make color palette for copy number state --------------------------------
cnv_state_colors <- c("gain" = "#e41a1c", "loss" = "#377eb8", "neutral" = "white")

# make colors for cell type -----------------------------------------------
## set color for different immune cell types
colors_immune_cell_types <- c("Macrophages M1" = "#ffff99", ## yellow
                              "Macrophages M1&M2" = "#fdbf6f", ## light orange
                              "Macrophages M2" = "#ff7f00", ## dark orange
                              "cDC1" = "#33a02c", ## dark green
                              "Myeloid lineage immune cells" = "#b2df8a", ## light green
                              "CD4+ T-cells" = "#1f78b4", ## dark blue
                              "CD8+ T-cells" = "#6a3d9a", ## dark purple
                              "NK cells" = "#cab2d6", ## light purple
                              "CD8+ T-cells & CD4+ T-cells" =  "#a6cee3", ## light blue
                              "B-cells" = "#fb9a99", ## pink
                              "Plasma cells" = "#e31a1c", ## red
                              "Unknown" = "grey50")

# color for cnv -----------------------------------------------------------
PuBu_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuBu")
PuRd_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuRd")

cna_state_colors <- c("Deep Loss" = PuBu_colors[9],
                      "Shallow Loss" = PuBu_colors[5],
                      "Neutral" = PuBu_colors[3],
                      "Low Gain" = PuRd_colors[5],
                      "High Gain" = PuRd_colors[9],
                      "Not Available" = "grey50")
copy_number_colors <-  c("0 Copies" = PuBu_colors[9],
                         "1 Copy" = PuBu_colors[5],
                         "2 Copies" = PuBu_colors[3],
                         "3 Copies" = PuRd_colors[5], 
                         "4 Copies" = PuRd_colors[7],
                         ">4 Copies" = PuRd_colors[9],
                         "Not Available" = "grey50")

save_pheatmap_pdf <- function(x, filename, width=6, height=6) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}



