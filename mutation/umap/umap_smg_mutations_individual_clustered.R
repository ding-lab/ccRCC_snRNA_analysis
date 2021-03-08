# Yige Wu @WashU Aug 2020
## 2020-09-03 removed the silent mutations

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
## set 10XMapping processing run id to input !0XMapping result later
mut_mapping_run_id <- "20200219.v1"
## input 10xmapping result
snRNA_mutation_df <- fread("./Resources/Analysis_Results/mutation/unite_10xmapping/20200303.v1/10XMapping.20200303.v1.tsv", data.table = F)
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input umap data
umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample/20200717.v1/Barcode2MetaData.20200717.v1.tsv")

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
for (snRNA_aliquot_id_tmp in unique(umap_df$aliquot)) {
  aliquot_show <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == snRNA_aliquot_id_tmp]
  file2write <- paste0(dir_out, aliquot_show, ".png")
  
  if (file.exists(file2write)) {
    next()
  }
  ## input barcodes with mapped varaint alleles and reference alleles
  mutation_map_tab.var <- snRNA_mutation_df %>%
    filter(aliquot == snRNA_aliquot_id_tmp) %>%
    filter(Is_Silent == F) %>%
    filter(gene_symbol %in% ccRCC_SMGs) %>%
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
  
  ### create read type, distinguish variant allele and reference allele
  plot_data_df$read_type <- "NA"
  plot_data_df$read_type[plot_data_df$barcode %in% mutation_map_tab.ref$barcode] <- "Ref"
  plot_data_df$read_type[plot_data_df$barcode %in% mutation_map_tab.var$barcode] <- "Var"
  
  ### create cell text, distinguish variant allele and reference allele
  plot_data_df <- plot_data_df %>%
    mutate(cell_text = paste0(str_split_fixed(string = mutation, pattern = "-Var", n = 2)[,1], "(", value, ")"))
  
  ### order the data frame so that cells mapped with variant allele will come on top
  plot_data_df <- rbind(plot_data_df[plot_data_df$read_type == "NA",],
                 plot_data_df[plot_data_df$read_type == "Ref",],
                 plot_data_df[plot_data_df$read_type == "Var",])
  
  ## make color palette for different read types
  colors_read_type <- c("#E31A1C", "#33A02C", "grey70")
  names(colors_read_type) <- c("Var", "Ref", "NA")
  ## ggplot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type != "Var",], mapping = aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 0.5, size = 0.3)
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type == "Var",], mapping = aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 0.8, size = 1)
  p <- p + scale_color_manual(values = colors_read_type)
  p <- p + ggtitle(paste0("Mapping of Mutations in SMG genes in ",  aliquot_show))
  # p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
  p <- p + geom_text_repel(data = plot_data_df[plot_data_df$cell_text != "(NA)",], mapping = aes(UMAP_1, UMAP_2, label = cell_text))
  p <- p +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "top")
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  p
  png(file2write, width = 800, height = 900, res = 150)
  print(p)
  dev.off()
}

