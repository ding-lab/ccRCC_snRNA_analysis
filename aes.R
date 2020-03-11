# especially for plotting #

# source ------------------------------------------------------------------
packages = c(
  "ggplot2",
  "RColorBrewer",
  "rcartocolor",
  "Polychrome"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# make color palette for variant class ------------------------------------
rcartocolor::display_carto_all()
cartocolors_df <- rcartocolor::cartocolors
variant_class_colors_green <- cartocolors_df[cartocolors_df$Name == "ag_GrnYl", "n4"][[1]]
names(variant_class_colors_green) <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", 'Splice_Site')
variant_class_colors_sunset <- cartocolors_df[cartocolors_df$Name == "ag_Sunset", "n3"][[1]]  
names(variant_class_colors_sunset) <- c("Missense_Mutation", "In_Frame_Ins", "Silent")
variant_class_colors <- c(variant_class_colors_green, variant_class_colors_sunset)
variant_class_colors

cartocolors_temps <- cartocolors_df[cartocolors_df$Name == "Temps", "n7"][[1]]
cartocolors_tropic <- cartocolors_df[cartocolors_df$Name == "Tropic", "n7"][[1]]

variant_class_colors <- c(cartocolors_temps[1:4], cartocolors_tropic[4], cartocolors_temps[c(5,7)])
names(variant_class_colors) <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", 'Splice_Site', "Silent", "Missense_Mutation", "In_Frame_Ins")
