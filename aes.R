# especially for plotting #

# source ------------------------------------------------------------------
packages = c(
  "ggplot2",
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer",
  "rcartocolor",
  "Polychrome"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}


# make color palette for each case ----------------------------------------
## input id meta data
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_cluster/write_path_to_seurat_objects_on_box/20200219.v1/Seurat_Object_Paths.20200219.v1.tsv", data.table = F)
uniq_case_ids <- unique(srat_paths$Case)
uniq_case_ids
### get unique color for each case
uniq_case_colors <- Polychrome::dark.colors(n = length(uniq_case_ids))
names(uniq_case_colors) <- uniq_case_ids


# make color palette for variant class ------------------------------------
# rcartocolor::display_carto_all()
cartocolors_df <- rcartocolor::cartocolors
cartocolors_temps <- cartocolors_df[cartocolors_df$Name == "Temps", "n7"][[1]]
cartocolors_tropic <- cartocolors_df[cartocolors_df$Name == "Tropic", "n7"][[1]]
variant_class_colors <- c(cartocolors_temps[1:4], cartocolors_tropic[4], cartocolors_temps[c(5,7)])
names(variant_class_colors) <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", 'Splice_Site', "Silent", "Missense_Mutation", "In_Frame_Ins")
