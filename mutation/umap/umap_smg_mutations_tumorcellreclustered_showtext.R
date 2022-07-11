# Yige Wu @WashU Aug 2020
## 2020-09-03 removed the silent mutations
## 2020-09-09 make the non-SMG mutations more transparent

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
## set 10XMapping processing run id to input !0XMapping result later
mut_mapping_run_id <- "20200219.v1"
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input 10xmapping result
snRNA_mutation_df <- fread("./Resources/Analysis_Results/mutation/unite_10xmapping/20200303.v1/10XMapping.20200303.v1.tsv", data.table = F)
## input umap data
umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_individual_tumorcellreclustered_on_katmai/20201127.v1/MetaData_TumorCellOnlyReclustered.20201127.v1.tsv")

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
snRNA_mutation_df$easy_id <- mapvalues(x = snRNA_mutation_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))

# remove the sample without mutation mapping result -----------------------
umap_df <- umap_df %>%
  filter(aliquot != "CPT0000890002")
aliquot2process <- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available & idmetadata_df$Sample_Type == "Tumor"]

# plot  by sample --------------------------------------------------
for (snRNA_aliquot_id_tmp in aliquot2process) {
  aliquot_show <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == snRNA_aliquot_id_tmp]
  
  ## input barcodes with mapped varaint alleles and reference alleles
  mutation_map_tab.var <- snRNA_mutation_df %>%
    filter(aliquot == snRNA_aliquot_id_tmp) %>%
    filter(Is_Silent == F) %>%
    filter(allele_type == "Var")
  
  if (snRNA_aliquot_id_tmp == "CPT0001260013") {
    mutation_map_tab.var <- mutation_map_tab.var %>%
      filter(gene_symbol != "FTH1")
  }
  
  mutation_map_tab.ref <- snRNA_mutation_df %>%
    filter(aliquot == snRNA_aliquot_id_tmp) %>%
    filter(Is_Silent == F) %>%
    filter(gene_symbol %in% mutation_map_tab.var$gene_symbol) %>%
    filter(allele_type == "Ref")
  
  ## make data frame for plotting
  plot_data_df <- umap_df %>%
    filter(orig.ident == snRNA_aliquot_id_tmp) %>%
    rename(barcode = barcode_tumorcellreclustered)
  
  ## merge with variant read info
  plot_data_df <- merge(plot_data_df, mutation_map_tab.var, by = c("barcode"), all.x = T)
  
  ### create read type, distinguish variant allele and reference allele
  plot_data_df$read_type <- "others"
  plot_data_df$read_type[plot_data_df$barcode %in% mutation_map_tab.var$barcode] <- "cells with the variant read(s)"
  table(plot_data_df$read_type)
  
  ### order the data frame so that cells mapped with variant allele will come on top
  plot_data_df <- rbind(plot_data_df[plot_data_df$read_type == "others",],
                        plot_data_df[plot_data_df$read_type == "cells with the variant read(s)",])
  
  ## add a column for drive genes
  plot_data_df <- plot_data_df %>%
    mutate(celltype_plot = ifelse(is.na(gene_symbol), "others",
                                  ifelse(gene_symbol %in% ccRCC_SMGs, "cells with the SMG mutation(s)", "others"))) %>%
    mutate(label_plot = ifelse(gene_symbol %in% ccRCC_SMGs, gene_symbol, NA))
  ## make color palette for different read types
  colors_celltype <- c("#E31A1C", "grey70")
  names(colors_celltype) <- c("cells with the SMG mutation(s)", "others")
  
  ## ggplot
  p <- ggplot(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = celltype_plot, label = label_plot))
  p <- p + geom_point_rast(alpha = 0.8, size = 0.5)
  p <- p + geom_text_repel(max.overlaps = Inf, min.segment.length = 0, box.padding = 0.5, size = 8, segment.alpha = 0.5, segment.size = 0.2)
  p <- p + scale_color_manual(values = colors_celltype)
  p <- p + ggtitle(paste0(aliquot_show))
  p <- p + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                 axis.text.y=element_blank(),axis.ticks=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 title = element_text(size = 20))
  p
  # file2write <- paste0(dir_out, aliquot_show, ".mut.png")
  # png(file2write, width = 800, height = 900, res = 150)
  # print(p)
  # dev.off()
  
  file2write <- paste0(dir_out, aliquot_show, ".mut.pdf")
  pdf(file2write, width = 4, height = 4, useDingbats = F)
  print(p)
  dev.off()
}

# count which samples have the most SMG mutation mapped -------------------
count_var_df <- snRNA_mutation_df %>%
  filter(allele_type == "Var") %>%
  filter(!Is_Silent) %>%
  filter(gene_symbol %in% ccRCC_SMGs) %>%
  group_by(easy_id) %>%
  summarise(count_var = n())


