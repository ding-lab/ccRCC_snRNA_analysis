# Yige Wu @WashU Oct 2020
## plot InferCNV subcluster mode outputs onto UMAP for individual sample (all cells)

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggrastr)

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input UMAP info per barcode
integrated_umap_df <- fread(input = "./Resources/Analysis_Results/fetch_data/fetchdata_34_ccRCC_samples_merged_katmai/20211005.v1/ccRCC.34Sample.Merged.Metadata.20211005.v1.tsv", data.table = F)
## input CNV state for known CNV genes
cnv_state_bycell_bygene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/annotate_barcode_with_cnv/annotate_barcode_with_gene_level_cnv_using_cnv_genes/20220131.v1/CNV_State_By_Gene_By_Barcode.20220131.v1.tsv")
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
genes2process <- knowncnvgenes_df$Gene_Symbol

# Loop: for each aliquot, input seurat object and infercnv output, plot important genes on UMAP ---------
gene_tmp <- "VHL"
# for (gene_tmp in "VHL") {
for (gene_tmp in genes2process) {
    
  ## extract current gene related results from the infercnv result data frame
  cnv_state_filtered_df <- cnv_state_bycell_bygene_df %>%
    filter(gene_symbol == gene_tmp)
  
  ## add CNV state to the barcode - UMAP coordidate data frame
  plotdata_df <- integrated_umap_df %>%
    mutate(barcode_individual = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1])
  plotdata_df <- merge(plotdata_df, cnv_state_filtered_df, 
                       by.x = c("orig.ident", "barcode_individual"), 
                       by.y = c("id_aliquot", "barcode_individual"),
                       all.x = T)
  unique(plotdata_df$orig.ident[is.na(plotdata_df$cna_value)])
  ## map CNV state value to text
  plotdata_df$cnv_cat <- map_infercnv_state2category(copy_state = plotdata_df$cna_value)
  plotdata_df$cnv_cat %>% table()
  
  ## make cells with CNV appear on top
  plotdata_df <- plotdata_df %>%
    arrange(factor(cnv_cat, levels = c("Not Available", "2 Copies", "1 Copy",  "3 Copies", ">4 Copies")))
  
  p <- ggplot() +
    geom_point_rast(data = plotdata_df, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat), alpha = 1, size = 0.3) +
    scale_color_manual(values = copy_number_colors)
  p <- p + ggtitle(paste0(gene_tmp, " Copy Number Status"))
  p <- p + theme_bw()
  p <- p + theme(panel.border = element_blank(), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "top")
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  
  file2write <- paste0(dir_out, gene_tmp, "_CNV_UMAP.png")
  png(file2write, width = 800, height = 900, res = 150)
  print(p)
  dev.off()
}

