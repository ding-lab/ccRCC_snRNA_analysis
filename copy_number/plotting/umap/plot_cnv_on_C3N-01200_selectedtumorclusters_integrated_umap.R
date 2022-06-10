# Yige Wu @WashU Mar 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "ggrastr",
  "ggplot2",
  "readxl"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode2UMAP
barcode2umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/integrate_C3N-01200_tumorcells/integrate_C3N-01200_tumorcells_selected_clusters/20220608.v1/C3N-01200.Tumorcells.Integrated.UMAP_data.20220608.v1.tsv")
## input known cnv genes
knowncnvgenes_df <- read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")

# pre-process -------------------------------------------------------------
## set genes to plot
genes2plot <- knowncnvgenes_df$Gene_Symbol
## set infercnv result direcotory
dir_infercnv_output <- "./Resources/snRNA_Processed_Data/InferCNV/outputs/"
path_infercnv_df <- data.frame(easy_id = paste0("C3N-01200-", c("T1", "T2", "T3")),
                               aliquot = c("CPT0075130004", "CPT0075140002", "CPT0075120002"),
                               dir_infercnv = paste0(dir_infercnv_output, c("Individual.20200305.v1/CPT0075130004/",
                                                                            "Individual.20200305.v1/CPT0075140002/",
                                                                            "Individual.20200305.v1/CPT0075120002/")))
## set colors
PuBu_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuBu")
PuRd_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuRd")
copy_number_colors <-  c("0 Copy" = PuBu_colors[9],
                         "1 Copy" = PuBu_colors[5],
                         "2 Copies" = "grey50",
                         "3 Copies" = PuRd_colors[5], 
                         "4 Copies" = PuRd_colors[7],
                         ">4 Copies" = PuRd_colors[9],
                         "Not Available" = "black")
copy_number_colors <- c("loss" = PuBu_colors[5],
                        "gain" = PuRd_colors[5],
                        "neutral" = "grey50")

# input CNV data ----------------------------------------------------
cnv_state_all_df <- NULL
for (i in 1:nrow(path_infercnv_df)) {
  ## input infercnv CNV state results
  obs_cnv_state_mat <- fread(input = paste0(path_infercnv_df$dir_infercnv[i], "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  ref_cnv_state_mat <- fread(input = paste0(path_infercnv_df$dir_infercnv[i], "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  cnv_state_df <- rbind(melt(obs_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))
  cnv_state_df <- cnv_state_df %>%
    filter(V1 %in% genes2plot) %>%
    mutate(easy_id = path_infercnv_df$easy_id[i]) %>%
    mutate(aliquot = path_infercnv_df$aliquot[i])
  cnv_state_all_df <- rbind(cnv_state_all_df, cnv_state_df)
}
rm(obs_cnv_state_mat)
rm(ref_cnv_state_mat)
table(cnv_state_all_df$aliquot)

# plot by gene ------------------------------------------------------------
gene_tmp <- "VHL"
for (gene_tmp in c("VHL", "SQSTM1")) {
  plot_data_df <- merge(x = barcode2umap_df %>%
                          mutate(barcode_individual = str_split_fixed(string = barcode_merged, pattern = "_", n = 2)[,1]), 
                        y = cnv_state_all_df %>%
                          filter(V1 == gene_tmp), by.x = c("orig.ident", "barcode_individual"), by.y = c("aliquot", "variable"), all.x = T)
  ## map CNV state value to text
  plot_data_df$cnv_cat <- map_infercnv_state2category(copy_state = plot_data_df$value)
  plot_data_df <- plot_data_df %>%
    mutate(cnv_cat_simple = ifelse(cnv_cat %in% c("0 Copy", "1 Copy"), "loss",
                                   ifelse(cnv_cat %in% c("3 Copies", "4 Copies", ">4 Copies"), "gain", "neutral")))
  
  p <- ggplot() +
    geom_point_rast(data = plot_data_df, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat_simple), shape = 16, alpha = 0.8, size = 0.2) +
    scale_color_manual(values = copy_number_colors)
  p <- p + theme_void()
  p <- p + theme(panel.border = element_blank(), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "none")
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  file2write <- paste0(dir_out, "C3N-01200.", gene_tmp, ".pdf")
  pdf(file2write, width = 2, height = 2, useDingbats = F)
  print(p)
  dev.off()
}

p <- ggplot() +
  geom_point_rast(data = plot_data_df, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat_simple), shape = 16, alpha = 0.8, size = 0.2) +
  scale_color_manual(values = copy_number_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = "Copy number status", nrow = 1))
p <- p + theme_void()
p <- p + theme(panel.border = element_blank(), 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
p <- p + theme(legend.position = "top")
p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())
file2write <- paste0(dir_out, "legend.", ".pdf")
pdf(file2write, width = 6, height = 4, useDingbats = F)
print(p)
dev.off()
