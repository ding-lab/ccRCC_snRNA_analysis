# Yige Wu @WashU Feb 2020
## for plotting the marker genes for integrated object, showing cell of origin

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
## load libraries
library(ggrepel)
# input dependencies ------------------------------------------------------
## input the paths for individual seurat object
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input 10xmapping result
snRNA_mutation_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/mutation/unite_10xmapping/20200303.v1/10XMapping.20200303.v1.tsv", data.table = F)

# plot  by sample --------------------------------------------------
for (snRNA_aliquot_id_tmp in unique(srat_paths$Aliquot)) {
  file2write <- paste(dir_out, snRNA_aliquot_id_tmp, ".SMG_Mutation_Mapping.", run_id, ".png", sep="")
  
  if (file.exists(file2write)) {
    next()
  }
  ## input barcodes with mapped varaint alleles and reference alleles
  mutation_map_tab.var <- snRNA_mutation_df %>%
    filter(aliquot == snRNA_aliquot_id_tmp) %>%
    filter(gene_symbol %in% ccRCC_SMGs) %>%
    filter(allele_type == "Var")
  
  mutation_map_tab.ref <- snRNA_mutation_df %>%
    filter(aliquot == snRNA_aliquot_id_tmp) %>%
    filter(gene_symbol %in% ccRCC_SMGs) %>%
    filter(allele_type == "Ref")
  
  ## input seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  ## get the coordinates for each cluster label
  p <- DimPlot(seurat_obj, reduction = "umap", label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  
  umap_tab <- FetchData(seurat_obj, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  umap_tab$barcode <- rownames(umap_tab)
  
  ## make data frame for plotting
  plot_data_df <- umap_tab
  
  ## merge with variant read info
  plot_data_df <- merge(plot_data_df, mutation_map_tab.var, by = c("barcode"), all.x = T)
  
  ### create read type, distinguish variant allele and reference allele
  plot_data_df$read_type <- "NA"
  plot_data_df$read_type[plot_data_df$barcode %in% mutation_map_tab.ref$barcode] <- "Ref"
  plot_data_df$read_type[plot_data_df$barcode %in% mutation_map_tab.var$barcode] <- "Var"
  
  ### create cell text, distinguish variant allele and reference allele
  plot_data_df <- plot_data_df %>%
    mutate(cell_text = str_split_fixed(string = mutation, pattern = "-Var", n = 2)[,1])
  
  ### order the data frame so that cells mapped with variant allele will come on top
  plot_data_df <- rbind(plot_data_df[plot_data_df$read_type == "NA",],
                 plot_data_df[plot_data_df$read_type == "Ref",],
                 plot_data_df[plot_data_df$read_type == "Var",])
  
  ## make color palette for different read types
  colors_read_type <- c("#E31A1C", "#33A02C", "grey70")
  names(colors_read_type) <- c("Var", "Ref", "NA")
  
  ## get the case id for this aliquot to show in the title
  case_id_tmp <- srat_paths$Case[srat_paths$Aliquot == snRNA_aliquot_id_tmp]
  
  ## ggplot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type != "Var",], mapping = aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 0.5, size = 0.3)
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type == "Var",], mapping = aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 0.8, size = 1)
  p <- p + scale_color_manual(values = colors_read_type)
  p <- p + ggtitle(paste0("Case: ", case_id_tmp, "   Aliquot: ",  snRNA_aliquot_id_tmp), 
                   subtitle = paste0("Mapping of All Mutations"))
  p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
  p <- p + geom_text_repel(data = plot_data_df[plot_data_df$cell_text != "",], mapping = aes(UMAP_1, UMAP_2, label = cell_text))
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
