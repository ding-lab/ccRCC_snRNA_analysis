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
library(ggridges)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input CNV state for known CNV genes
cnv_state_bycell_bygene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/annotate_barcode_with_cnv/annotate_barcode_with_gene_level_cnv_using_cnv_genes/20220202.v1/CNV_State_By_Gene_By_Barcode.20220202.v1.tsv")
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
genes2process <- knowncnvgenes_df$Gene_Symbol

# Loop: for each aliquot, input seurat object and infercnv output, plot important genes on UMAP ---------
gene_tmp <- "VHL"
for (gene_tmp in "VHL") {
# for (gene_tmp in genes2process) {
    
  ## extract current gene related results from the infercnv result data frame
  cnv_state_filtered_df <- cnv_state_bycell_bygene_df %>%
    filter(gene_symbol == gene_tmp) %>%
    filter(id_aliquot != "CPT0002270013")
  
  ## add CNV state to the barcode - UMAP coordidate data frame
  plotdata_df <- merge(cnv_state_filtered_df, 
                       barcode2celltype_df %>%
                         select(orig.ident, individual_barcode, Cell_group5), 
                       by.x = c("id_aliquot", "barcode_individual"), 
                       by.y = c("orig.ident", "individual_barcode"),
                       all.x = T)
  plotdata_df <- plotdata_df %>%
    filter(Cell_group5 != "Unknown")
  
  p <- ggplot()
  p <- p + geom_density_ridges(data = plotdata_df, mapping = aes(x = cna_value, y = Cell_group5, fill = Cell_group5), rel_min_height = 0.01, scale = 0.9, color = NA)
  p <- p + scale_fill_manual(values = colors_cellgroup5)
  p <- p + ggtitle(paste0(gene_tmp, " Copy Number Ratio"))
  p <- p + xlim(c(quantile(plotdata_df$cna_value, 0.01)-0.75, quantile(plotdata_df$cna_value, 0.99)+.75))
  p <- p + xlab("copy number ratio")
  p <- p + theme_bw()
  p <- p + theme(axis.title.y = element_blank(), axis.text = element_text(color = "black"),
                 legend.position = "none")
  # p <- p + guides(fill = guide_legend(nrow = 3, title = NULL))
  
  file2write <- paste0(dir_out, gene_tmp, ".predHMMi6.ridgeline.png")
  png(file2write, width = 500, height = 400, res = 150)
  print(p)
  dev.off()
  
  # p <- ggplot()
  # p <- p + geom_density_ridges(data = plotdata_df, mapping = aes(x = cna_value, y = Cell_group5, fill = Cell_group5), 
  #                              stat = "binline", bins = 20,
  #                              rel_min_height = 0.01, scale = 0.9)
  # p <- p + ggtitle(paste0(gene_tmp, " Copy Number Ratio"))
  # p <- p + xlim(c(quantile(plotdata_df$cna_value, 0.01)-0.5, quantile(plotdata_df$cna_value, 0.99)+0.5))
  # p <- p + xlab("copy number ratio")
  # p <- p + theme_bw()
  # p <- p + theme(axis.title.y = element_blank(), axis.text = element_text(color = "black"))
  # 
  # file2write <- paste0(dir_out, gene_tmp, ".predHMMi6.ridgeline2.png")
  # png(file2write, width = 800, height = 400, res = 150)
  # print(p)
  # dev.off()
}

