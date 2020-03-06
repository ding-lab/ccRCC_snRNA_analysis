# Yige Wu @WashU Feb 2020
## getting the expressin of CNv-related DEGs in cells with or without CNVs per sample

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the paths for individual seurat object
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input DEG united list
deg_by_chr_region_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/findallmarkers/unite_degs_by_chr_region/20200228.v1/FindMarkers.Wilcox.ExpectedCNV_vs_Neutral.20200228.v1.tsv", data.table = F)
## input the average expression data 
plot_data_sup_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/findallmarkers/get_exp_of_deg_by_chr_region/20200228.v1/CNV_Marker_Average_Exp_Data.20200228.v1.tsv", data.table = F)
## input the genetic downstream table
genetic_alt_downstream_genes <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_ccrcc_genetic_event_downstream_genes/20200227.v1/ccRCC_Genetic_Event_Downstream_Genes.20200227.v1.tsv", data.table = F)
## set chr_region to process
chr_region_tmp <- "3p"

# make the data frame for plotting ----------------------------------------
## filter by chromosomal region
plot_data_df <- plot_data_sup_df %>%
  filter(chr_region == chr_region_tmp) %>%
  filter(id != "Other")

# ggplot ------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plot_data_df, mapping = aes(x = features.plot, y = aliquot, size = pct.exp, fill = avg.exp.scaled), shape = 21)
p <- p + scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                               high = "red", space = "Lab" )
p <- p + facet_grid(id~., scales = "free", space = "free", drop = T)
p <- p + theme(axis.text.x = element_text(angle = 90, size = 8))
p
png(filename = paste0(dir_out, chr_region_tmp, "_Recurrent_DEG_Expression.", run_id, ".png"), width = 800, height = 300)
print(p)
dev.off



