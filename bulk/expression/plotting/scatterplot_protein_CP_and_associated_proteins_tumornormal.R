# Yige Wu @WashU March 2020
## reference: https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggpubr)
library(ggrastr)

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input protein data
protein_raw_df <- fread("./Resources/Bulk_Processed_Data/Protein/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)
## input bulk meta data for the entire set
meta_df <- fread("./Resources/Bulk_Processed_Data/Meta_Data/cptac-metadata.csv", data.table = F)
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# specify parameters ------------------------------------------------------
## input filtered deg
gene_x <- "CP"
# genes_y <- c("HIF1A", "MXI1", "SREBF2", "VEGFA")
genes_y <- c( "VEGFA")
genes_plot <- c(gene_x, genes_y)
## input id variables
id_colnames <- c("Index", "Proteins", "ReferenceIntensity")

# preprocess --------------------------------------------------------------
# get the aliquot IDs for bulk --------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma")
aliquotids_plot <- c(metadata_filtered_df$Specimen.Label.tumor, metadata_filtered_df$Specimen.Label.normal)
aliquotids_plot <- aliquotids_plot[!is.na(aliquotids_plot)]
# aliquotids_plot <- meta_df$Specimen.Label[!is.na(meta_df$Type)]
## make the matrix to plot the heatmap
protein_raw_filtered_df <- protein_raw_df %>%
  filter(Index %in% genes_plot)
protein_raw_filtered_mat <- protein_raw_filtered_df %>%
  select(aliquotids_plot)
rownames(protein_raw_filtered_mat) <- protein_raw_filtered_df$Index
## deduct reference intensity
protein_plot_mat <- as.matrix(protein_raw_filtered_mat) - as.vector(protein_raw_filtered_df$ReferenceIntensity)
protein_plot_df <- cbind(protein_raw_filtered_df[, id_colnames], protein_plot_mat)
plot_data_dfx <- melt(data = protein_plot_df[protein_plot_df$Index == gene_x,], id.vars = id_colnames) 
plot_data_dfx$Sample_Type <- mapvalues(x = plot_data_dfx$variable, from = meta_df$Specimen.Label, to = as.vector(meta_df$Type))

for (gene_y in genes_y) {
  # make plot data ----------------------------------------------------------
  plot_data_dfy <- melt(data = protein_plot_df[protein_plot_df$Index == gene_y,], id.vars = id_colnames) 
  plot_data_df <- merge(x = plot_data_dfx, y = plot_data_dfy, by = c("variable"), suffixes = c(".x", ".y"))
  
  # make scatterplot --------------------------------------------------------
  p <- ggplot()
  # p <- p + geom_point(data = plot_data_df, mapping = aes(x = x, y = y))
  p <- ggscatter(data = plot_data_df, x = "value.x", y = "value.y", color = "Sample_Type", size = 0.5,
                 add = "reg.line",  # Add regressin line
                 # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE # Add confidence interval
  )
  # Add correlation coefficient
  p <- p + stat_cor(aes(color = Sample_Type),
                    method = "pearson",
                    label.x = quantile(x = plot_data_df$value.x, probs = 0.001, na.rm = T))
  p <- p + xlab(paste0(gene_x, " protein abundance"))
  p <- p + ylab(paste0(gene_y, " protein abundance"))
  p <- p + theme_classic(base_size = 12)
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14, color = "black"))
  p

  # save scatterplot --------------------------------------------------------
  # file2write <- paste0(dir_out, gene_x, "~", gene_y,".png")
  # png(file2write, width = 800, height = 800, res = 150)
  # print(p)
  # dev.off()

  file2write <- paste0(dir_out, gene_x, "~", gene_y,".pdf")
  pdf(file2write, width = 3, height = 3, useDingbats = F)
  print(p)
  dev.off()
  
  # # make scatterplot --------------------------------------------------------
  # p <- ggplot()
  # # p <- p + geom_point(data = plot_data_df, mapping = aes(x = x, y = y))
  # p <- ggscatter(data = plot_data_df, x = "value.x", y = "value.y", color = "Sample_Type", size = 0.5,
  #                add = "reg.line",  # Add regressin line
  #                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  #                conf.int = TRUE # Add confidence interval
  # )
  # # Add correlation coefficient
  # p <- p + stat_cor(method = "pearson", 
  #                   label.x = quantile(x = plot_data_df$value.x, probs = 0.001, na.rm = T), 
  #                   label.y = quantile(x = plot_data_df$value.y, probs = 0.99, na.rm = T))
  # p <- p + xlab(paste0(gene_x, " protein abundance"))
  # p <- p + ylab(paste0(gene_y, " protein abundance"))
  # p <- p + theme_classic(base_size = 9)
  # p <- p + theme(legend.position = "top", title = element_text(size = 9))
  # p <- p + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8, color = "black"))
  # p
  # 
  # # save scatterplot --------------------------------------------------------
  # # file2write <- paste0(dir_out, gene_x, "~", gene_y,".png")
  # # png(file2write, width = 800, height = 800, res = 150)
  # # print(p)
  # # dev.off()
  # 
  # file2write <- paste0(dir_out, gene_x, "~", gene_y,".pdf")
  # pdf(file2write, width = 2, height = 2.3, useDingbats = F)
  # print(p)
  # dev.off()
  
  # # make scatterplot --------------------------------------------------------
  # p <- ggplot()
  # p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = value.x, y = value.y))
  # p <- ggscatter(data = plot_data_df, x = "value.x", y = "value.y", color = "Sample_Type", size = 0.5,
  #                add = "reg.line",  # Add regressin line
  #                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  #                conf.int = TRUE # Add confidence interval
  # )
  # # Add correlation coefficient
  # p <- p + stat_cor(method = "pearson", 
  #                   label.x = quantile(x = plot_data_df$value.x, probs = 0.001, na.rm = T), 
  #                   label.y = quantile(x = plot_data_df$value.y, probs = 0.99, na.rm = T))
  # p <- p + xlab(paste0(gene_x, " protein abundance"))
  # p <- p + ylab(paste0(gene_y, " protein abundance"))
  # p <- p + theme_classic(base_size = 9)
  # p <- p + theme(legend.position = "top", title = element_text(size = 9))
  # p <- p + theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8, color = "black"))
  # p
  # 
  # # save scatterplot --------------------------------------------------------
  # # file2write <- paste0(dir_out, gene_x, "~", gene_y,".png")
  # # png(file2write, width = 800, height = 800, res = 150)
  # # print(p)
  # # dev.off()
  # 
  # file2write <- paste0(dir_out, gene_x, "~", gene_y,".pdf")
  # pdf(file2write, width = 2, height = 2.3, useDingbats = F)
  # print(p)
  # dev.off()
  
  
}
