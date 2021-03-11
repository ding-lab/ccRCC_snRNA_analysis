# Yige Wu @WashU Mar 2021

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
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210305.v1/meta_data.20210305.v1.tsv")
## input 10xmapping result
snRNA_mutation_df <- fread("./Resources/Analysis_Results/mutation/unite_10xmapping_snatac_mapping_heatmap_outputs/20201224.v1/10XMapping.snATAC.mapping_heatmap.20201224.v1.tsv", data.table = F)

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

# plot  by sample --------------------------------------------------
easyids_process <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$snATAC_used]
for (easyid_tmp in easyids_process) {
# for (easyid_tmp in c("C3L-00416-T2")) {
    
  snRNA_aliquot_id_tmp <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU == easyid_tmp]
  path_umap <- paste0("../ccRCC_snATAC/Resources/snATAC_Processed_Data/Barcode_Annotation/UMAP/Individual_Sample/", snRNA_aliquot_id_tmp, "_UMAP_data.tsv")
  if (!file.exists(path_umap)) {
    next()
  }
  umap_df <- fread(data.table = F, input = path_umap)
  
  ## input barcodes with mapped varaint alleles and reference alleles
  mutation_map_tab.var <- snRNA_mutation_df %>%
    filter(aliquot == snRNA_aliquot_id_tmp) %>%
    filter(Is_Silent == F) %>%
    filter(allele_type == "Var")

  mutation_map_tab.ref <- snRNA_mutation_df %>%
    filter(aliquot == snRNA_aliquot_id_tmp) %>%
    filter(Is_Silent == F) %>%
    filter(gene_symbol %in% mutation_map_tab.var$gene_symbol) %>%
    filter(allele_type == "Ref")
  
  ## make data frame for plotting
  plot_data_df <- umap_df %>%
    rename(barcode = Barcode)
  
  ## merge with variant read info
  plot_data_df <- merge(plot_data_df, mutation_map_tab.var, by = c("barcode"), all.x = T)
  
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
  p <- p + ggtitle(paste0("Mapping of Mutations in ",  easyid_tmp, " snATAC-seq data"))
  p <- p + geom_text_repel(data = plot_data_df[!is.na(plot_data_df$gene_symbol) & plot_data_df$Driver_Gene_Mutation == "FALSE",],
                           mapping = aes(UMAP_1, UMAP_2, label = gene_symbol, colour = Driver_Gene_Mutation, size = Driver_Gene_Mutation), segment.alpha = 0.3)
  p <- p + geom_text_repel(data = plot_data_df[!is.na(plot_data_df$gene_symbol) & plot_data_df$Driver_Gene_Mutation == "TRUE",],
                           mapping = aes(UMAP_1, UMAP_2, label = gene_symbol, colour = Driver_Gene_Mutation, size = Driver_Gene_Mutation), segment.alpha = 1)
  p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey30"))
  p <- p + scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3))
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
  file2write <- paste0(dir_out, easyid_tmp, ".mut.png")
  png(file2write, width = 800, height = 900, res = 150)
  print(p)
  dev.off()
  
}
