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
## input the barcode to CNV state info
barcode2cnv_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/annotate_barcode_cnv/20200228.v1/Expected_CNV_State_By_Chr_Region_By_Barcode.20200228.v1.tsv", data.table = F)
## input DEG united list
deg_by_chr_region_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/findallmarkers/unite_degs_by_chr_region/20200228.v1/FindMarkers.Wilcox.ExpectedCNV_vs_Neutral.20200228.v1.tsv", data.table = F)
## input the average expression data 
plot_data_sup_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/findallmarkers/get_exp_of_deg_by_chr_region/20200228.v1/CNV_Marker_Average_Exp_Data.20200228.v1.tsv", data.table = F)

# plot by chr region ------------------------------------------------------

for (chr_region_tmp in unique(deg_count_by_chr_region_df$chr_region)) {
  ## count the frequency of each gene associated with each chr_region
  deg_count_by_chr_region_df <- deg_by_chr_region_df %>%
    select(gene, chr_region) %>%
    filter(chr_region == chr_region_tmp) %>%
    table() %>%
    as.data.frame() %>%
    filter(Freq > 0) %>%
    arrange(desc(Freq))
  
  ## filter by chromosomal region
  plot_data_df <- plot_data_sup_df %>%
    filter(chr_region == chr_region_tmp) %>%
    filter(id != "Other")
  
  ## order the genes
  plot_data_df$features.plot <- factor(plot_data_df$features.plot, levels = as.vector(deg_count_by_chr_region_df$gene))
  ## ggplot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = features.plot, y = aliquot, size = pct.exp, fill = avg.exp.scaled), shape = 21)
  p <- p + scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                                  high = "red", space = "Lab" )
  p <- p + facet_grid(id~., scales = "free", space = "free", drop = T)
  p <- p + theme(axis.text.x = element_text(angle = 90, size = 8))
  p
  png(filename = paste0(dir_out, chr_region_tmp, ".DEG_Expression.", run_id, ".png"), width = 2000, height = 300)
  stop("")
}