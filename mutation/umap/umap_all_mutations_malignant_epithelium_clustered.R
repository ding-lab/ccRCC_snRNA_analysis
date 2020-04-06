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

# input dependencies ------------------------------------------------------
## set 10XMapping processing run id to input !0XMapping result later
mut_mapping_run_id <- "20200219.v1"
## input seurat processing info to input individual seurat object later
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200207.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object

# plot  by sample --------------------------------------------------
for (snRNA_aliquot_id_tmp in unique(seurat_summary2process$Aliquot)) {
  file2write <- paste(dir_out, snRNA_aliquot_id_tmp, ".All_Mutation_Mapping.", run_id, ".png", sep="")
  
  if (file.exists(file2write)) {
    next()
  }
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  ## get the coordinates for each cluster label
  p <- DimPlot(seurat_obj, reduction = "umap", label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  
  umap_tab <- FetchData(seurat_obj, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  umap_tab$barcode <- rownames(umap_tab)
  
  ## input barcodes with mapped varaint alleles and reference alleles
  mutation_map_tab <- fread(input = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/10Xmapping/outputs/", 
                                           mut_mapping_run_id, "/",
                                           snRNA_aliquot_id_tmp, "/", 
                                           snRNA_aliquot_id_tmp, "_mapping_heatmap_0.txt"), data.table = F)
  mutation_map_tab <- mutation_map_tab %>%
    mutate(gene_symbol = str_split_fixed(string = Mutatation, pattern = "-", n = 3)[,1])
  mutation_map_tab.m <- melt(mutation_map_tab, id.vars = c("gene_symbol", "Mutatation"))
  mutation_map_tab.m <- mutation_map_tab.m %>%
    mutate(allele_type = str_split_fixed(string = Mutatation, pattern = "-", n = 3)[,3])
  
  mutation_map_tab.var <- mutation_map_tab.m %>%
    filter(allele_type == "Var") %>%
    filter(!is.na(value) & value > 0) %>%
    rename(barcode = variable)
  
  mutation_map_tab.ref <- mutation_map_tab.m %>%
    filter(allele_type == "Ref") %>%
    filter(!is.na(value) & value > 0) %>%
    rename(barcode = variable)
  
  ## make data frame for plotting
  plot_data_df <- umap_tab
  
  ### create read type, distinguish variant allele and reference allele
  plot_data_df$read_type <- "NA"
  plot_data_df$read_type[plot_data_df$barcode %in% mutation_map_tab.ref$barcode] <- "Ref"
  plot_data_df$read_type[plot_data_df$barcode %in% mutation_map_tab.var$barcode] <- "Var"
  table(plot_data_df$read_type)
  
  ### order the data frame so that cells mapped with variant allele will come on top
  plot_data_df <- rbind(plot_data_df[plot_data_df$read_type == "NA",],
                 plot_data_df[plot_data_df$read_type == "Ref",],
                 plot_data_df[plot_data_df$read_type == "Var",])
  
  ## make color palette for different read types
  colors_read_type <- c("#E31A1C", "#33A02C", "grey70")
  names(colors_read_type) <- c("Var", "Ref", "NA")
  
  ## get the case id for this aliquot to show in the title
  case_id_tmp <- seurat_summary2process$Case[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  
  ## ggplot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type != "Var",], mapping = aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 0.5, size = 0.3)
  p <- p + geom_point(data = plot_data_df[plot_data_df$read_type == "Var",], mapping = aes(x = UMAP_1, y = UMAP_2, color=read_type), alpha = 0.8, size = 0.8)
  p <- p + scale_color_manual(values = colors_read_type)
  p <- p + ggtitle(paste0("Case: ", case_id_tmp, "   Aliquot: ",  snRNA_aliquot_id_tmp), 
                   subtitle = paste0("Mapping of All Mutations"))
  p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = ident))
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
