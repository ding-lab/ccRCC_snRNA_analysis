# Yige Wu @WashU Aug 2020
## 2020-09-03 removed the silent mutations
## 2020-09-09 make the non-SMG mutations more transparent

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
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
## set 10XMapping processing run id to input !0XMapping result later
mut_mapping_run_id <- "20200219.v1"
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input 10xmapping result
snRNA_mutation_df <- fread("./Resources/Analysis_Results/mutation/unite_10xmapping/20200303.v1/10XMapping.20200303.v1.tsv", data.table = F)
## input umap data
umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetch_data_by_individual_sample/20200717.v1/Barcode2MetaData.20200717.v1.tsv")
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20200902.v1/scrublet.run20200902_adj_cutoff.united_outputs.tsv", data.table = F)

# remove the silent mutations ---------------------------------------------
snRNA_mutation_df <- snRNA_mutation_df %>%
  mutate(AA_Change = str_split_fixed(string = mutation, pattern = "-", n = 3)[,2]) %>%
  mutate(Ref_Allele = str_split_fixed(string = AA_Change, pattern = '[0-9]', n = 2)[,1])
snRNA_mutation_df$Alt_Allele <- sapply(X = snRNA_mutation_df$AA_Change, FUN = function(x) {
  text_vec <- str_split(string = x, pattern = '[0-9]')[[1]]
  return(text_vec[length(text_vec)])
})  
snRNA_mutation_df <- snRNA_mutation_df %>%
  mutate(Is_Silent = (Alt_Allele == Ref_Allele))

# remove the sample without mutation mapping result -----------------------
umap_df <- umap_df %>%
  filter(aliquot != "CPT0000890002")

# plot  by sample --------------------------------------------------
for (snRNA_aliquot_id_tmp in "CPT0001260013") {
  aliquot_show <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == snRNA_aliquot_id_tmp]
  
  ## input barcodes with mapped varaint alleles and reference alleles
  mutation_map_tab.var <- snRNA_mutation_df %>%
    filter(aliquot == snRNA_aliquot_id_tmp) %>%
    filter(Is_Silent == F) %>%
    filter(gene_symbol != "FTH1") %>%
    filter(allele_type == "Var")
  
  mutation_map_tab.ref <- snRNA_mutation_df %>%
    filter(aliquot == snRNA_aliquot_id_tmp) %>%
    filter(Is_Silent == F) %>%
    filter(gene_symbol %in% mutation_map_tab.var$gene_symbol) %>%
    filter(allele_type == "Ref")
  
  ## make data frame for plotting
  plot_data_df <- umap_df %>%
    filter(aliquot == snRNA_aliquot_id_tmp) %>%
    rename(barcode = individual_barcode)
  
  ## merge with variant read info
  plot_data_df <- merge(plot_data_df, mutation_map_tab.var, by = c("barcode"), all.x = T)
  plot_data_df <- merge(x = plot_data_df, y = barcode2scrublet_df, by.x = c("aliquot.x", "barcode"), by.y = c("Aliquot", "Barcodes"), all.x = T)
  
  ### create read type, distinguish variant allele and reference allele
  plot_data_df$read_type <- "NA"
  plot_data_df$read_type[plot_data_df$barcode %in% mutation_map_tab.var$barcode] <- "Var"
  table(plot_data_df$read_type)
  
  ### order the data frame so that cells mapped with variant allele will come on top
  plot_data_df <- rbind(plot_data_df[plot_data_df$read_type == "NA",],
                 plot_data_df[plot_data_df$read_type == "Var",])
  
  ## add a column for drive genes
  plot_data_df <- plot_data_df %>%
    mutate(mutation_cat = ifelse(gene_symbol %in% ccRCC_drivers, "SMG",
                                 ifelse(gene_symbol == "C8orf76", "Shared", "Others"))) %>%
    mutate(Driver_Gene_Mutation = ifelse(gene_symbol %in% ccRCC_drivers, "TRUE", "FALSE")) %>%
    mutate(read_type_text = ifelse(read_type == "NA", "others", "cells with the variant read(s)"))
  ## make color palette for different read types
  colors_read_type <- c("#E31A1C", "grey70")
  names(colors_read_type) <- c("cells with the variant read(s)", "others")
  
  ## ggplot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type == "NA",], mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.8, size = 0.5, color = colors_read_type["others"])
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type == "Var" & !(plot_data_df$gene_symbol %in% ccRCC_SMGs),], mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.5, size = 0.8, color = colors_read_type["cells with the variant read(s)"])
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type == "Var" & (plot_data_df$gene_symbol %in% ccRCC_SMGs),], mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.8, size = 1, color = colors_read_type["cells with the variant read(s)"])
  p <- p + ggtitle(paste0("Mapping of Mutations in ",  aliquot_show))
  p <- p + geom_text_repel(data = plot_data_df[!is.na(plot_data_df$gene_symbol),],
                           mapping = aes(UMAP_1, UMAP_2, label = gene_symbol, colour = Driver_Gene_Mutation, size = Driver_Gene_Mutation), max.overlaps = 30)
  p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey50"))
  p <- p + scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3))
  p <- p + facet_grid(cols = vars(predicted_doublet))
  p <- p +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "none")
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  p
  file2write <- paste0(dir_out, aliquot_show, ".mut.png")
  png(file2write, width = 2000, height = 900, res = 150)
  print(p)
  dev.off()
  
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type == "NA",], mapping = aes(x = UMAP_1, y = UMAP_2, color = read_type_text), alpha = 0.5, size = 0.3)
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type == "Var",], mapping = aes(x = UMAP_1, y = UMAP_2, color = read_type_text), alpha = 0.5, size = 0.8)
  p <- p + scale_color_manual(values = colors_read_type)
  p <- p + ggtitle(paste0("Mapping of Mutations in ",  aliquot_show))
  p <- p +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p <- p + labs(colour = "")
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme(legend.position = "top")
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  p
  file2write <- paste0(dir_out, aliquot_show, ".mut.legend.png")
  png(file2write, width = 800, height = 900, res = 150)
  print(p)
  dev.off()
  
  ## ggplot
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$read_type == "NA",], mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.8, size = 0.5, color = colors_read_type["others"])
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$read_type == "Var" & !(plot_data_df$gene_symbol %in% ccRCC_SMGs),], mapping = aes(x = UMAP_1, y = UMAP_2), 
                      alpha = 1, size = 2, color = colors_read_type["cells with the variant read(s)"])
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$read_type == "Var" & (plot_data_df$gene_symbol %in% ccRCC_SMGs),], mapping = aes(x = UMAP_1, y = UMAP_2), 
                      alpha = 1, size = 2.5, color = colors_read_type["cells with the variant read(s)"])
  # p <- p + geom_text_repel(data = plot_data_df[!is.na(plot_data_df$gene_symbol),], 
  #                          mapping = aes(UMAP_1, UMAP_2, label = gene_symbol, colour = mutation_cat, size = Driver_Gene_Mutation), fontface = "italic")
  p <- p + geom_text_repel(data = plot_data_df[!is.na(plot_data_df$gene_symbol) & plot_data_df$mutation_cat %in% c("SMG", "Shared"),], 
                           mapping = aes(UMAP_1, UMAP_2, label = gene_symbol, colour = mutation_cat, size = Driver_Gene_Mutation), fontface = "italic", force = 2)
  p <- p + scale_color_manual(values = c("SMG" = "black", "Shared" = "black", "FALSE" = "grey50"))
  p <- p + scale_size_manual(values = c("TRUE" = 10, "FALSE" = 7))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p <- p + theme(axis.title = element_text(size = 20))
  p <- p + theme(legend.position = "none")
  p
  file2write <- paste0(dir_out, aliquot_show, ".mut.nolgend.pdf")
  pdf(file2write, width = 8, height = 8, useDingbats = F)
  print(p)
  dev.off()
  
  p <- ggplot()
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$read_type == "NA",], mapping = aes(x = UMAP_1, y = UMAP_2, color = read_type_text), alpha = 0.5, size = 0.3)
  p <- p + geom_point_rast(data = plot_data_df[plot_data_df$read_type == "Var",], mapping = aes(x = UMAP_1, y = UMAP_2, color = read_type_text), alpha = 0.5, size = 0.8)
  p <- p + scale_color_manual(values = colors_read_type)
  p <- p +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p <- p + labs(colour = "")
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  p <- p + theme(legend.position = "bottom", legend.text = element_text(size = 20), legend.title = element_text(size = 25))
  p
  file2write <- paste0(dir_out, aliquot_show, ".mut.lgend.pdf")
  pdf(file2write, width = 8, height = 8, useDingbats = F)
  print(p)
  dev.off()
}
