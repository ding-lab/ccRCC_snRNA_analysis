# Yige Wu @WashU Nov 2019
## for plotting the marker genes for integrated object, showing cell of origin

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# set aliquot ids to be processed -----------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191024.v1.tsv")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191021.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds")) %>%
  mutate(Paht_deg_table = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                 "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                 "/", Aliquot, FACS, ".DEGs.Pos.txt"))
seurat_summary2process$Path_seurat_object


# plot dotplot by sample --------------------------------------------------
snRNA_aliquot_id_tmp <- "CPT0001220012"
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  umap_tab <- FetchData(seurat_obj, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  umap_tab$barcode_int <- rownames(umap_tab)
  umap_tab <- umap_tab %>%
    filter(orig.ident %in% snRNA_aliquot_id_tmp) %>%
    mutate(barcode = str_split_fixed(string = barcode_int, pattern = "_", n = 2)[,1])
  
  # input 10xmapping results ---------------------------------------------
  mutation_map_tab <- fread(input = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/10Xmapping/outputs/", snRNA_aliquot_id_tmp, "/", snRNA_aliquot_id_tmp, "_mapping_heatmap_0.txt"), data.table = F)
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
  
  genes2plot <- mutation_map_tab.var %>%
    select(gene_symbol) %>%
    unique() 
  genes2plot <- genes2plot$gene_symbol
  
  for (gene_tmp in genes2plot) {
    mutation_map_tab.var.tmp <- mutation_map_tab.var %>%
      filter(value > 0) %>%
      filter(gene_symbol == gene_tmp) %>%
      mutate(mutation = paste0(gene_tmp, "-", str_split_fixed(string = Mutatation, pattern = "-", n = 3)[,2]))
    
    mutation_map_tab.ref.tmp <- mutation_map_tab.ref %>%
      filter(value > 0) %>%
      filter(gene_symbol == gene_tmp)
    
    tab2p <- umap_tab
    tab2p$read_type <- mapvalues(x = tab2p$barcode, from = mutation_map_tab.var.tmp$barcode, to = mutation_map_tab.var.tmp$mutation)
    tab2p$read_type[tab2p$barcode %in% mutation_map_tab.ref.tmp$barcode] <- "Ref"
    tab2p$read_type[tab2p$read_type %in% tab2p$barcode] <- "NA"
    table(tab2p$read_type)
    
    tab2p <- rbind(tab2p[tab2p$read_type == "NA",],
                   tab2p[tab2p$read_type == "Ref",],
                   tab2p[tab2p$read_type %in% unique(mutation_map_tab.var.tmp$mutation),])
    
    colors_read_type <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1")[1:length(unique(mutation_map_tab.var.tmp$mutation))], "green", "grey70")
    names(colors_read_type) <- c(unique(mutation_map_tab.var.tmp$mutation), "Ref", "NA")
    
    # plot a UMAP plot for copy number metric ---------------------------------
    p <- ggplot(tab2p, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(color=read_type), alpha = 1, size = 0.3) +
      scale_color_manual(values = colors_read_type) +
      ggtitle(paste0(snRNA_aliquot_id_tmp, "_", gene_tmp, "_Mutation_Status")) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank())
    file2write <- paste(dir_out, snRNA_aliquot_id_tmp,".", gene_tmp, ".FeaturePlot_Mutation_Count.", run_id, ".png", sep="")
    png(file2write, width = 1000, height = 800, res = 150)
    print(p)
    dev.off()
  }
}
