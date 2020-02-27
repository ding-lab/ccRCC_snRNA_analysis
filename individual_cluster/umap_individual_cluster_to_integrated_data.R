# Yige Wu @WashU Feb 2020
## for plotting individual sample clusters onto the integrated object to confirm the cell group for the individual cluster

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
## input integrated data UMAP coordinates
integrated_umap_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/20200212.v3/30_aliquot_integration.20200212.v3.umap_data.tsv", data.table = F)
### format the barcode for the sake of merging
integrated_umap_df <- integrated_umap_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1])
### get label coordinates
label_data_df <- integrated_umap_df %>%
  group_by(ident) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# plot umap per aliquot ---------------------------------------------------
for (aliquot_tmp in unique(seurat_summary2process$Aliquot)) {
  ## input barcode to cluster info
  path_file <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/inputs/annotations_file/Individual.20200207.v1/", aliquot_tmp, ".Barcode_Annotation.txt")
  barcode2cluster_tmp <- fread(input = path_file, data.table = F)
  
  ## create sub output directory
  dir_out_sub1 <- paste0(dir_out, aliquot_tmp, "/")
  dir.create(dir_out_sub1)
  
  ## loop per cluster
  for (cluster_tmp in unique(barcode2cluster_tmp$V2)) {
    ### annotate the data frame for plotting to highlight the barcodes 
    plot_data_df <- integrated_umap_df
    ### get barcodes to highlight
    highlight_barcodes <- barcode2cluster_tmp$V1[barcode2cluster_tmp$V2 == cluster_tmp]
    ### note barcodes to highlight in the data frame
    plot_data_df$hlt_barcode <- (plot_data_df$barcode_individual %in% highlight_barcodes)
    ### make ggplot
    p <- ggplot()
    p <- p + geom_point(data = plot_data_df[!(plot_data_df$hlt_barcode),], 
                        mapping = aes(x = UMAP_1, y = UMAP_2), 
                        alpha = 0.3, size = 0.3, color = "grey50")
    p <- p + geom_point(data = plot_data_df[(plot_data_df$hlt_barcode),], 
                        mapping = aes(x = UMAP_1, y = UMAP_2), 
                        alpha = 0.8, size = 0.8, color = "red")
    p <- p + ggtitle(paste0("Aliquot: ",  aliquot_tmp, " Cluster: ", cluster_tmp), 
                     subtitle = paste0("Mapping onto integrated dataset"))
    p <- p + geom_text_repel(data = label_data_df, mapping = aes(UMAP_1, UMAP_2, label = ident), size = 5, colour = "white")
    # p <- p + theme_bw()
    p <- p + ggplot2::theme(legend.position = "right",
                            legend.text = element_text(size = 9))
    p <- p + ggplot2::theme(panel.border = element_blank(), 
                            panel.background = element_rect(fill = "black", colour = "black"),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank())
    p <- p + ggplot2::theme(axis.line=element_blank(),
                            axis.text.x=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank())
    p
    ### save ggplot
    file2write <- paste0(dir_out_sub1, "aliquot_", aliquot_tmp, ".cluster_", cluster_tmp, ".mapped_onto_integrated_dataset.png")
    png(file2write, width = 2000, height = 2100, res = 150)
    print(p)
    dev.off()
  }
}

