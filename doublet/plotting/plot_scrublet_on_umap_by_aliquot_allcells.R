# Yige Wu @WashU Set 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type per barcode table
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20200902.v1/scrublet.run20200902_adj_cutoff.united_outputs.tsv", data.table = F)
## input UMAP info per barcode
umap_df <- fread(input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample//20200406.v1/individual_cluster_meta_data.20200406.v1.tsv", data.table = F)
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# plot for cell group----------------------------------------------------------
for (id_aliquot_tmp in unique(umap_df$aliquot)) {
  aliquot_show <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == id_aliquot_tmp]
  ## filter down to one sample
  plotdata_df <- umap_df %>%
    filter(aliquot == id_aliquot_tmp)
  barcode2scrublet_aliquot_df <- barcode2scrublet_df %>%
    filter(Aliquot == id_aliquot_tmp)
  plotdata_df <- merge(plotdata_df,
                       barcode2scrublet_aliquot_df,
                       by.x = c("individual_barcode"), by.y = c("Barcodes"), all.x = T)
  
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, 
                      mapping = aes(x = UMAP_1, y = UMAP_2, color = predicted_doublet),
                      alpha = 1, size = 0.05)
  p <- p + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50"))
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p <- p + ggtitle(label = paste0(aliquot_show, " Cell Type"))
  p <- p + theme(legend.position = "top")
  p
  file2write <- paste0(dir_out, aliquot_show, ".png")
  png(filename = file2write, width = 800, height = 900, res = 150)
  print(p)
  dev.off()
  # file2write <- paste0(dir_out, aliquot_show, ".pdf")
  # pdf(file2write, width = 8, height = 9, useDingbats = F)
  # print(p)
  # dev.off()
}
