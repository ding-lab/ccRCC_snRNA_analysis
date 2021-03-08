# especially for plotting #

# source ------------------------------------------------------------------
packages = c(
  "ggplot2",
  "ggrepel",
  "ComplexHeatmap",
  "ggrastr",
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
colors_cellgroup5 <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")[c(4, 1, 3, 6)]
colors_cellgroup5 <- c(colors_cellgroup5, "grey50")
names(colors_cellgroup5) <- c("Tumor cells", "Normal epithelial cells", "Immune",  "Stroma",
                              "Unknown")
# swatch(colors_cellgroup5)

# make color palette for 13 cell groups -----------------------------------
# colors_cellgroup13 <- c(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[c(4, 1, 3)], ggthemes::colorblind_pal()(8)[c(1, 2, 3, 5, 6, 7, 8)], RColorBrewer::brewer.pal(n = 12, name = "Paired")[c(3, 5)], "grey50")
colors_cellgroup13 <- c("#E7298A", "#1B9E77", "#7570B3","#000000", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00","#CC79A7", "#B2DF8A", "#FB9A99","grey50")
names(colors_cellgroup13) <- c("Tumor cells", "Normal epithelial cells", "Immune others", "B-cells", "CD4+ T-cells", "CD8+ T-cells", "Macrophages", "DC", "NK cells","Endothelial cells", "Myofibroblasts", "Fibroblasts","Unknown")
# swatch(colors_cellgroup13)

# make color palette for 14 cell groups with transitional cells-----------------------------------
colors_cellgroup14 <- c(colors_cellgroup13, RColorBrewer::brewer.pal(n = 12, name = "Paired")[c(12)])
names(colors_cellgroup14) <- c(names(colors_cellgroup13), "EMT tumor cells")
colors_cellgroup14

# make color palette for variant class ------------------------------------
# rcartocolor::display_carto_all()
cartocolors_temps <- cartocolors_df[cartocolors_df$Name == "Temps", "n7"][[1]]
cartocolors_tropic <- cartocolors_df[cartocolors_df$Name == "Tropic", "n7"][[1]]
variant_class_colors <- c(cartocolors_temps[1:4], 
                          cartocolors_tropic[4], 
                          cartocolors_temps[c(5,6,7)],
                          "white")
names(variant_class_colors) <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", 'Splice_Site', 
                                 "Silent", 
                                 "Missense_Mutation", "In_Frame_Ins", "In_Frame_Del",
                                 "None")
colors_variant_class_sim <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[c(1,2,3,4,5)]
names(colors_variant_class_sim) <- c("Missense", "Truncation", "In_Frame_Ins", "Missense,Truncation", "In_Frame_Del")

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



